use proc_macro::TokenStream;


#[proc_macro_attribute]
pub fn metis_func(_input: TokenStream, annotated_item: TokenStream) -> TokenStream {

    let mut signature = annotated_item.clone().into_iter();
    let mut found = false;
    let mut idx = 0;
    for token in &mut signature {
        idx += 1;
        if &token.to_string() == "fn" {
            found = true;
            break;
        }
    }
    assert_eq!(found, true, "must be attached to an fn");

    if idx > 30 {
        eprintln!("Warning, metis_func searching further than expected");
    }

    let fn_name = signature.next().unwrap().to_string();

    let link_func = format!("libmetis__{}", fn_name);
    let mut output: TokenStream = quote::quote!{
        #[export_name = #link_func]
        #[allow(non_snake_case)]
    }.into();

    output.extend(annotated_item);

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

            let args_dec: Vec<_> = item_fn.sig.inputs.iter().map(|i| match i {
                    syn::FnArg::Typed(syn::PatType{pat, ..}) => pat,
                    _ => panic!("metis extern functions can't take self")
                }
            ).collect();

            let args = &item_fn.sig.inputs;
            let ret_type = &item_fn.sig.output;

            ext_signature.sig.ident = syn::Ident::new(&link_func_str, proc_macro2::Span::call_site()); 
            ext_signature.vis = syn::Visibility::Inherited;
            let extern_decl: TokenStream = quote::quote! {
                extern "C" {
                    #ext_signature
                }
                #[inline(always)]
                pub unsafe fn #rust_func ( #args ) #ret_type {
                    #link_func ( #(#args_dec),*)
                }
            }.into();
            output.extend(extern_decl)
        } else {
            panic!("all metis declarations must be functions");
        }
    }

    // eprintln!("{}", output);
    // panic!();

    output
}
