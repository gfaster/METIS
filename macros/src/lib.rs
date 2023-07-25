use proc_macro::TokenStream;
use quote::format_ident;
use syn::{ItemFn, Token};

#[proc_macro_attribute]
pub fn metis_func(_input: TokenStream, annotated_item: TokenStream) -> TokenStream {
    let mut item_fn = syn::parse_macro_input!(annotated_item as ItemFn);
    let fn_name = &item_fn.sig.ident;

    let mut foreign: syn::ForeignItemFn = syn::ForeignItemFn {
        attrs: item_fn.attrs.clone(),
        vis: syn::parse_quote!(pub(crate)),
        sig: item_fn.sig.clone(),
        semi_token: Token![;](proc_macro2::Span::call_site()),
    };
    foreign.sig.ident = format_ident!("libmetis__{}", fn_name);
    foreign.sig.ident.set_span(fn_name.span());
    foreign.sig.abi = None;
    foreign.sig.unsafety = None;

    let link_func = format!("libmetis__{}", fn_name);

    item_fn.sig.unsafety = Some(Token!(unsafe)(proc_macro2::Span::call_site()));

    // foreign block is so ab_tests can still call out to those functions
    let output: TokenStream = quote::quote! {
        #[cfg(feature = "dual_link")]
        extern "C" {
            #foreign
        }

        #[cfg_attr(not(feature = "dual_link"), export_name = #link_func)]
        #[allow(non_snake_case)]
        #item_fn
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

            ext_signature.sig.ident =
                syn::Ident::new(&link_func_str, proc_macro2::Span::call_site());
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
            #[cfg_attr(not(feature = "dual_link"), ignore)]
            #[test]
            #original_fn_decl {
                use super::#c_target as #target_name;
                #(#body_c)*
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

    item.attrs.push(syn::parse_quote!(#[test]));
    item.sig.output = syn::ReturnType::Default;
    item.block.stmts = syn::parse_quote! {
        #original_fn_decl {
            let #target_name = #c_target;
            #(#body_c)*
        }
        #rust_fn_decl {
            let #target_name = #rust_target;
            #(#body_rust)*
        }

        assert_eq!(rust(), original());
    };

    let output = quote::quote! {
        #item
    };

    TokenStream::from(output)
    // eprintln!("{}", output);
    // panic!();
}
