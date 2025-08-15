mod common;

mod gpmetis {
    use super::common::*;

    #[test]
    fn basic() {
        let p = Params::new(GpMetis, TestGraph::Elt4, 10);
        p.call_assert_eq();
    }

    #[test]
    #[ignore = "todo: fixme"]
    fn basic_all() {
        Params::new(GpMetis, "", 10).call_suite_assert_eq();
    }

    #[test]
    fn c_preserved() {
        let mut p = Params {
            ptype: Some(Ptype::Pmetis),
            overrides: "*:c",
            ..Params::new(GpMetis, "", 10)
        };
        p.call_suite_assert_eq();
        p.ptype = Some(Ptype::Kmetis);
        p.call_suite_assert_eq();
    }

    #[test]
    #[ignore = "todo: fixme"]
    fn basic_rb() {
        Params {
            ptype: Some(Ptype::Pmetis),
            overrides: "*:c",
            ..Params::new(GpMetis, TestGraph::Elt4, 10)
        }.call_suite_assert_eq();
    }

    #[test]
    fn help_msg() {
        Params {
            help: true,
            ..Params::new(GpMetis, TestGraph::Elt4, 0)
        }.call_assert_eq();
    }

    #[test]
    #[ignore = "broken in original"]
    fn minconn_vol() {
        Params {
            minconn: true,
            objtype: Some(Objtype::Vol),
            ..Params::new(GpMetis, "", 10)
        }.call_suite_assert_eq();
    }
}
