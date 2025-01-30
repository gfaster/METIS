#![allow(clippy::let_and_return)]

use proc_macro::TokenStream;
use proc_macro2::{Ident, Span};
use quote::format_ident;
use syn::{ItemFn, Token};

#[proc_macro_attribute]
pub fn metis_func(input: TokenStream, annotated_item: TokenStream) -> TokenStream {
    let item_fn = syn::parse_macro_input!(annotated_item as ItemFn);
    let mut pfx = "libmetis__";
    if !input.is_empty() {
        let directive = syn::parse_macro_input!(input as Ident);
        let directive_text = directive.to_string();
        if directive_text != "no_pfx" {
            panic!("unknown directive {directive_text} (try \"no_pfx\")");
        }
        pfx = "";
    }
    metis_func_normal(item_fn, pfx)
}

fn metis_func_normal(mut impl_fn: ItemFn, pfx: &str) -> TokenStream {
    let fn_name = &impl_fn.sig.ident;

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


    let mut dispatch_sig = impl_fn.sig.clone();
    dispatch_sig.ident = format_ident!("{fn_name}");
    dispatch_sig.ident.set_span(fn_name.span());
    dispatch_sig.unsafety = Some(Token!(unsafe)(Span::call_site()));
    dispatch_sig.abi = Some(syn::Abi {
        extern_token: Token![extern](Span::call_site()),
        name: Some(syn::LitStr::new("C", Span::call_site())),
    });
    let dispatch_vis = impl_fn.vis.clone();
    let mut dispatch_rs_call = format_ident!("{rs_link_func}");
    dispatch_rs_call.set_span(impl_fn.sig.ident.span());
    let dispatch_c_call = foreign.sig.ident.clone();
    let dispatch_c_sym_lit = {
        let sym = format!("{dispatch_c_call}\0");
        let sym = std::ffi::CStr::from_bytes_with_nul(sym.as_bytes()).expect("c symbol has invalid name");
        proc_macro2::Literal::c_string(sym)
    };

    let dispatch_args_decl: Vec<_> = impl_fn
        .sig
        .inputs
        .iter()
        .map(|i| match i {
            syn::FnArg::Typed(syn::PatType { pat, .. }) => pat,
            _ => panic!("metis extern functions can't take self"),
        }).cloned()
        .collect();

    let dispatch_args = impl_fn.sig.inputs.clone();

    // let underscore = Token![_](proc_macro2::Span::call_site());
    // let underscores = vec![underscore; dispatch_args_decl.len()];
    let mut dispatch_return = impl_fn.sig.output.clone();
    if matches!(dispatch_return, syn::ReturnType::Default) {
        dispatch_return = syn::ReturnType::Type(
            Token![->](Span::call_site()),
            syn::Type::Tuple(syn::TypeTuple {
                paren_token: syn::token::Paren::default(),
                elems: syn::punctuated::Punctuated::new(),
            }).into()
        )
    }
    // let dispatch_args_decl = quote::quote! { #(#dispatch_args_decl),* };

    let dispatch_fn_ptr = quote::quote! {
        // fn () -> ()
        unsafe extern "C" fn(#dispatch_args) #dispatch_return
    };

    // We do resolve functions so that we can link to them from the original source
    let resolve_name = format_ident!("resolve_{fn_name}");

    let dispatch = quote::quote! {
        #[export_name = #cannonical_link_func]
        #[allow(non_snake_case)]
        #dispatch_vis #dispatch_sig {
            let actual: #dispatch_fn_ptr = #resolve_name();
            unsafe { actual(#(#dispatch_args_decl),*) }
            // type FnPtr = fn ();
        }

        #[doc(hidden)]
        #[no_mangle]
        pub extern "C" fn #resolve_name() -> #dispatch_fn_ptr {
            static #dispatch_lib_ident: crate::dyncall::ICall = unsafe { crate::dyncall::ICall::new(
                #dispatch_c_sym_lit,
                &crate::dyncall::LIBMETIS,
                {
                    let f: #dispatch_fn_ptr = #dispatch_rs_call;
                    std::mem::transmute(f)
                }
            )};
            unsafe { std::mem::transmute(#dispatch_lib_ident.get()) }
        }
    };

    let prev_span = impl_fn.sig.ident.span();
    impl_fn.sig.ident = format_ident!("{rs_link_func}");
    impl_fn.sig.ident.set_span(prev_span);
    impl_fn.sig.abi = Some(syn::Abi { 
        extern_token: Token![extern](Span::call_site()), 
        name: Some(syn::LitStr::new("C", Span::call_site())) 
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
