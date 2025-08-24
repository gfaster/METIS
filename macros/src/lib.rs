#![allow(clippy::let_and_return)]

use proc_macro::TokenStream;
use proc_macro2::{Ident, Span};
use quote::{format_ident, ToTokens};
use syn::{punctuated::Punctuated, spanned::Spanned, ItemFn, Signature, Token};

struct ErrorCollector {
    inner: TokenStream,
}

impl ErrorCollector {
    #[cold]
    fn push_error(&mut self, err: syn::Error) {
        self.inner.extend(TokenStream::from(err.to_compile_error()));
    }

    #[cold]
    fn push<T: std::fmt::Display>(&mut self, span: Span, msg: T) {
        self.push_error(syn::Error::new(span, msg));
    }

    fn push_res(&mut self, res: Result<(), syn::Error>) {
        if let Err(e) = res {
            self.push_error(e);
        }
    }

    fn finish(mut self, rem: impl Into<TokenStream>) -> TokenStream {
        if self.inner.is_empty() {
            return rem.into()
        }

        self.inner.extend(rem.into());

        self.inner
    }
}

/// Function that is implemented in both C and Rust
///
/// Accepted arguments are:
///
/// ### `no_pfx`
/// exported name is not prefixed with `libmetis__`. This is used for unmangled names, i.e. api
/// names that start with `METIS_`.
///
/// ### `no_c`
/// Never use the C implementation. If it doesn't exist or exists and is ifunc'd, then it will
/// always call the Rust version. Does nothing if the C impl exists but is not ifunc'd. 
///
/// The C prototype of course must exist if it is to be called by C code. Remember to add the entry
/// in `rename.h` as well
#[proc_macro_attribute]
pub fn metis_func(input: TokenStream, annotated_item: TokenStream) -> TokenStream {
    let mut pfx = "libmetis__";
    let mut no_c = false;
    let mut errors = ErrorCollector { inner: TokenStream::new() };
    if !input.is_empty() {
        let directives = syn::parse_macro_input!(input with Punctuated<Ident, Token![,]>::parse_terminated);

        for directive in directives {
            let text = directive.to_string();
            match &*text {
                "no_pfx" => pfx = "",
                "no_c" => no_c = true,
                _ => {
                    errors.push_error(
                    syn::Error::new(directive.span(),
                        format!("unknown directive `{text}` (valid directives: `no_pfx` or `no_c`)"))
                    );
                }
            }
        }
    }

    if !errors.inner.is_empty() {
        return errors.finish(annotated_item)
    }

    let item_fn = syn::parse_macro_input!(annotated_item as ItemFn);

    if no_c {
        metis_func_no_c_impl(item_fn, pfx, errors)
    } else {
        metis_func_normal(item_fn, pfx, errors)
    }
}

fn ensure_extern_c(sig: &Signature, msg: &'static str) -> Result<(), syn::Error> {
    let Some(abi) = sig.abi.as_ref() else {
        return Err(syn::Error::new(sig.span(), msg))
    };
    let Some(lit) = abi.name.as_ref() else {
        return Err(syn::Error::new(sig.span(), msg))
    };
    if lit.value() != "C" {
        return Err(syn::Error::new(sig.span(), msg))
    }
    Ok(())
}

fn ensure_unsafe(sig: &Signature, msg: &'static str) -> Result<(), syn::Error> {
    let Some(_unsafety) = sig.unsafety.as_ref() else {
        return Err(syn::Error::new(sig.span(), msg))
    };
    Ok(())
}

fn metis_func_no_c_impl(impl_fn: ItemFn, pfx: &str, mut errors: ErrorCollector) -> TokenStream {
    errors.push_res(ensure_extern_c(&impl_fn.sig, "no_c metis_func must be `extern \"C\"`"));
    errors.push_res(ensure_unsafe(&impl_fn.sig, "metis_func must be `unsafe`"));
    let export_name = format!("{pfx}{}", impl_fn.sig.ident);
    // let resolve_name_exported = format_ident!("resolve_{fn_name}");
    let dispatch_fn_ptr = fn_pointer_of_sig(&impl_fn.sig, &mut errors);
    let impl_ident = &impl_fn.sig.ident;
    let resolve_name_exported = format_ident!("resolve_{fn_name}", fn_name = impl_fn.sig.ident);


    let res = quote::quote! {
        #[export_name = #export_name]
        #[allow(dead_code, non_snake_case)]
        #impl_fn

        #[doc(hidden)]
        #[no_mangle]
        pub extern "C" fn #resolve_name_exported() -> #dispatch_fn_ptr {
            #impl_ident
        }
    };
    errors.finish(res)
}

fn dispatch_sig(base_sig: &Signature) -> Signature {
    let mut dispatch_sig = base_sig.clone();
    dispatch_sig.ident = base_sig.ident.clone();
    dispatch_sig.unsafety = Some(Token!(unsafe)(Span::call_site()));
    // always need extern C because in tests, dispatch is what the ifunc resolves to
    dispatch_sig.abi = Some(syn::Abi {
        extern_token: Token![extern](Span::call_site()),
        name: Some(syn::LitStr::new("C", Span::call_site())),
    });
    dispatch_sig
}

fn fn_pointer_of_sig(sig: &Signature, errors: &mut ErrorCollector) -> proc_macro2::TokenStream {
    if !sig.generics.params.is_empty() {
        for ty in sig.generics.type_params() {
            errors.push_error(syn::Error::new(ty.span(), "cannot create function pointer out of function with type params"));
        }
        for cnst in sig.generics.const_params() {
            errors.push_error(syn::Error::new(cnst.span(), "cannot create function pointer out of function with const params"));
        }
        for lt in sig.generics.lifetimes() {
            errors.push_error(syn::Error::new(lt.span(), "cannot create function pointer out of function with explicit lifetime params"));
        }
    }
    if let Some(var) = &sig.variadic {
        errors.push_error(syn::Error::new(var.span(), "cannot create function pointer out of variadic function"));
    }

    let mut bare_arg = |arg: &syn::FnArg| {
        match arg {
            syn::FnArg::Receiver(receiver) => {
                errors.push(receiver.span(), "function pointer cannot have reciever");
                syn::BareFnArg {
                    attrs: receiver.attrs.clone(),
                    name: None,
                    ty: (*receiver.ty).clone(),
                }
            },
            syn::FnArg::Typed(pat_type) => {
                syn::BareFnArg {
                    attrs: pat_type.attrs.clone(),
                    name: None,
                    ty: (*pat_type.ty).clone(),
                }
            },
        }
    };

    let inputs = sig.inputs.pairs().map(|p| match p {
        syn::punctuated::Pair::Punctuated(arg, comma) => {
            syn::punctuated::Pair::Punctuated(bare_arg(arg), *comma)
        },
        syn::punctuated::Pair::End(arg) => {
            syn::punctuated::Pair::End(bare_arg(arg))
        },
    }).collect();

    syn::TypeBareFn {
        lifetimes: None,
        unsafety: sig.unsafety,
        abi: sig.abi.clone(),
        fn_token: sig.fn_token,
        paren_token: sig.paren_token,
        inputs,
        variadic: None,
        output: sig.output.clone(),
    }.into_token_stream()
}

fn metis_func_normal(mut impl_fn: ItemFn, pfx: &str, mut errors: ErrorCollector) -> TokenStream {
    let fn_name = &impl_fn.sig.ident;
    let attrs = impl_fn.attrs.clone();

    let mut foreign: syn::ForeignItemFn = syn::ForeignItemFn {
        attrs: impl_fn.attrs.clone(),
        vis: syn::parse_quote!(pub),
        sig: impl_fn.sig.clone(),
        semi_token: Token![;](Span::call_site()),
    };
    foreign.sig.ident = format_ident!("c__{pfx}{}", fn_name);
    foreign.sig.ident.set_span(fn_name.span());
    foreign.sig.abi = None;
    foreign.sig.unsafety = None;

    let cannonical_link_func = format!("{pfx}{fn_name}");

    let rs_link_func = format!("rs__{pfx}{}", fn_name);
    let dispatch_lib_ident = format_ident!("__LIBRARY_DISPATCH_{fn_name}");

    let dispatch_sig = dispatch_sig(&impl_fn.sig);
    let dispatch_sig_ident = &dispatch_sig.ident;
    let dispatch_vis = impl_fn.vis.clone();
    let mut dispatch_rs_call = format_ident!("{rs_link_func}");
    dispatch_rs_call.set_span(impl_fn.sig.ident.span());
    let dispatch_c_call = foreign.sig.ident.clone();
    let dispatch_c_sym_lit = {
        let sym = format!("{dispatch_c_call}\0");
        let sym =
            std::ffi::CStr::from_bytes_with_nul(sym.as_bytes()).expect("c symbol has invalid name");
        proc_macro2::Literal::c_string(sym)
    };

    let dispatch_fn_ptr = fn_pointer_of_sig(&dispatch_sig, &mut errors);

    let dispatch_args_decl: Vec<_> = impl_fn
        .sig
        .inputs
        .iter()
        .map(|i| match i {
            syn::FnArg::Typed(syn::PatType { pat, .. }) => pat,
            _ => panic!("metis extern functions can't take self"),
        })
        .cloned()
        .collect();

    // let underscore = Token![_](proc_macro2::Span::call_site());
    // let underscores = vec![underscore; dispatch_args_decl.len()];
    // let dispatch_args_decl = quote::quote! { #(#dispatch_args_decl),* };

    // let dispatch_fn_ptr = quote::quote! {
    //     // fn () -> ()
    //     unsafe extern "C" fn(#dispatch_args) #dispatch_return
    // };

    // We do resolve functions so that we can link to them from the original source
    let resolve_name = format_ident!("resolve_static_{fn_name}");
    let resolve_name_exported = format_ident!("resolve_{fn_name}");

    let dispatch = quote::quote! {
        #[export_name = #cannonical_link_func]
        #[allow(non_snake_case)]
        #(#attrs)*
        #dispatch_vis #dispatch_sig {
            let actual: #dispatch_fn_ptr = #resolve_name();
            unsafe { actual(#(#dispatch_args_decl),*) }
        }

        #[doc(hidden)]
        #[no_mangle]
        pub extern "C" fn #resolve_name_exported() -> #dispatch_fn_ptr {
            if cfg!(test) {
                #dispatch_sig_ident
            } else {
                #resolve_name()
            }
        }

        #[doc(hidden)]
        fn #resolve_name() -> #dispatch_fn_ptr {
            #![deny(unused)]

            static #dispatch_lib_ident: crate::dyncall::ICall = unsafe { crate::dyncall::ICall::new(
                #dispatch_c_sym_lit,
                &crate::dyncall::LIBMETIS,
                {
                    let f: #dispatch_fn_ptr = #dispatch_rs_call;
                    std::mem::transmute(f)
                }
            )};
            #[cfg(test)]
            let actual: std::ptr::NonNull<std::ffi::c_void> = {
                thread_local!{
                    static EPOCH: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
                    static CURRENT: std::cell::Cell<std::ptr::NonNull<std::ffi::c_void>> =
                        const { std::cell::Cell::new(std::ptr::NonNull::dangling()) };
                }
                let actual_epoch = crate::dyncall::ICall::epoch();
                if actual_epoch != EPOCH.get() {
                    let p = #dispatch_lib_ident.get_dynamic();
                    CURRENT.set(p);
                    EPOCH.set(actual_epoch);
                    p
                } else {
                    CURRENT.get()
                }
            };

            #[cfg(not(test))]
            let actual: std::ptr::NonNull<std::ffi::c_void> = {
                #dispatch_lib_ident.get()
            };

            unsafe { std::mem::transmute::<std::ptr::NonNull<std::ffi::c_void>, #dispatch_fn_ptr>(actual) }
        }
    };

    let prev_span = impl_fn.sig.ident.span();
    impl_fn.sig.ident = format_ident!("{rs_link_func}");
    impl_fn.sig.ident.set_span(prev_span);
    impl_fn.sig.abi = Some(syn::Abi {
        extern_token: Token![extern](Span::call_site()),
        name: Some(syn::LitStr::new("C", Span::call_site())),
    });

    impl_fn.sig.unsafety = Some(Token!(unsafe)(Span::call_site()));

    // foreign block is so ab_tests can still call out to those functions
    let output: TokenStream = quote::quote! {
        // #[cfg(any(feature = "dual_link", feature = "no_rs"))]
        // extern "C" {
        //     #[allow(clippy::too_many_arguments)]
        //     #foreign
        // }

        #dispatch

        // #[export_name = #rs_link_func]
        #[allow(non_snake_case, clippy::too_many_arguments)]
        #[doc(hidden)]
        #impl_fn
    }
    .into();

    // eprintln!("{}", output);

    output
}

#[proc_macro_attribute]
pub fn metis_decl(_input: TokenStream, annotated_item: TokenStream) -> TokenStream {
    let input = syn::parse_macro_input!(annotated_item as syn::ItemForeignMod);

    let mut output = TokenStream::new();

    for item in input.items {
        if let syn::ForeignItem::Fn(item_fn) = item {
            let mut ext_signature = item_fn.clone();
            let link_func = quote::format_ident!("libmetis__{}", item_fn.sig.ident);
            let link_func_str = format!("libmetis__{}", item_fn.sig.ident);
            let rust_func = item_fn.sig.ident;

            let args_dec: Vec<_> = item_fn
                .sig
                .inputs
                .iter()
                .map(|i| match i {
                    syn::FnArg::Typed(syn::PatType { pat, .. }) => pat,
                    _ => panic!("metis extern functions can't take self"),
                })
                .collect();

            let args = &item_fn.sig.inputs;
            let ret_type = &item_fn.sig.output;

            ext_signature.sig.ident = syn::Ident::new(&link_func_str, rust_func.span());
            ext_signature.vis = syn::Visibility::Inherited;
            let extern_decl: TokenStream = quote::quote! {
                extern "C" {
                    #[allow(non_snake_case)]
                    #ext_signature
                }
                #[inline(always)]
                pub unsafe fn #rust_func ( #args ) #ret_type {
                    #link_func ( #(#args_dec),*)
                }
            }
            .into();
            output.extend(extern_decl)
        } else {
            panic!("all metis declarations must be functions");
        }
    }

    // eprintln!("{}", output);
    // panic!();

    output
}

#[proc_macro_attribute]
pub fn ab_test_basic(input: TokenStream, annotated_item: TokenStream) -> TokenStream {
    let target = syn::parse_macro_input!(input as syn::Path);
    let target_name = &target
        .segments
        .last()
        .expect("ab_test argument must be a valid path")
        .ident;

    let item = syn::parse_macro_input!(annotated_item as syn::ItemFn);
    let mod_name = item.sig.ident.clone();
    let fn_span = item.sig.ident.span();
    let sig = item.sig;

    let rust_target = target.clone();

    let c_target_name = format_ident!("libmetis__{}", target_name);
    let mut c_target = target.clone();
    c_target
        .segments
        .last_mut()
        .expect("ab_test argument must be a valid path")
        .ident = c_target_name;

    let body_rust = item.block.stmts.clone().into_iter();
    let body_c = item.block.stmts.into_iter();

    let mut original_fn_decl = sig.clone();
    original_fn_decl.ident = syn::Ident::new("original", fn_span);

    let mut rust_fn_decl = sig.clone();
    rust_fn_decl.ident = syn::Ident::new("rust", fn_span);

    let output = quote::quote! {
        mod #mod_name {
            #[cfg_attr(not(feature = "dual_link"), ignore = "requires dual_link")]
            #[test]
            #original_fn_decl {
                #[cfg(feature = "dual_link")]
                {
                    use super::#c_target as #target_name;
                    #(#body_c)*
                }
            }
            #[test]
            #rust_fn_decl {
                use super::#rust_target as #target_name;
                #(#body_rust)*
            }
        }
    }
    .into();

    // eprintln!("{}", output);
    // panic!();

    output
}

#[proc_macro_attribute]
pub fn ab_test_eq(input: TokenStream, annotated_item: TokenStream) -> TokenStream {
    let target = syn::parse_macro_input!(input as syn::Path);
    let target_name = &target
        .segments
        .last()
        .expect("ab_test argument must be a valid path")
        .ident;

    let mut item = syn::parse_macro_input!(annotated_item as syn::ItemFn);
    let fn_span = item.sig.ident.span();
    let rust_target = target.clone();

    let c_target_name = format_ident!("libmetis__{}", target_name);
    let mut c_target = target.clone();
    c_target
        .segments
        .last_mut()
        .expect("ab_test argument must be a valid path")
        .ident = c_target_name;

    let body_rust = item.block.stmts.clone().into_iter();
    let body_c = item.block.stmts.into_iter();

    let mut original_fn_decl = item.sig.clone();
    original_fn_decl.ident = syn::Ident::new("original", fn_span);
    let mut rust_fn_decl = item.sig.clone();
    rust_fn_decl.ident = syn::Ident::new("rust", fn_span);

    item.attrs
        .push(syn::parse_quote!(#[cfg_attr(not(feature = "dual_link"), ignore)]));
    item.attrs.push(syn::parse_quote!(#[test]));
    item.sig.output = syn::ReturnType::Default;
    item.block.stmts = syn::parse_quote! {
        #[cfg(feature = "dual_link")]
        {
            #original_fn_decl {
                #[allow(non_snake_case)]
                let #target_name = #c_target;
                #(#body_c)*
            }
            #rust_fn_decl {
                #[allow(non_snake_case)]
                let #target_name = #rust_target;
                #(#body_rust)*
            }

            assert_eq!(rust(), original());
        }
    };

    let output = quote::quote! {
        #item
    };

    TokenStream::from(output)
    // eprintln!("{}", output);
    // panic!();
}
