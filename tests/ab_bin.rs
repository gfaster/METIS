mod common;

mod gpmetis {
    use super::common::*;

    #[test]
    fn basic() {
        let p = Params::new(GpMetis, TestGraph::Elt4, 10);
        p.call_assert_eq();
    }

    #[test]
    fn basic_all() {
        Params::new(GpMetis, "", 10).call_suite_assert_eq();
    }

    #[test]
    fn c_preserved() {
        let mut p = Params {
            overrides: "*:c",
            ..Params::new(GpMetis, "", 10)
        };
        p.call_suite_assert_eq();
    }

    #[test]
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
    fn minconn_cut() {
        // I thought this was broken in the original??? but for some reason it works fine??
        Params {
            minconn: true,
            objtype: Some(Objtype::Cut),
            ..Params::new(GpMetis, "", 10)
        }.call_suite_assert_eq();
    }

    #[test]
    fn minconn_vol() {
        // If metis-normalized is built with assertions, this will fail. I expect that it has other
        // failures as well that are a bit more common.
        Params {
            minconn: true,
            objtype: Some(Objtype::Vol),
            ..Params::new(GpMetis, "", 10)
        }.call_suite_assert_eq();
    }

    #[test]
    #[ignore = "slow"]
    fn little_bit_of_everthing_pmetis() {
        let mut p = Params::new(GpMetis, "", 0);
        p.nooutput = true;
        // p.dbglvl = Some(metis::METIS_DBG_IPART | metis::METIS_DBG_MOVEINFO | metis::METIS_DBG_REFINE | metis::METIS_DBG_COARSEN | metis::METIS_DBG_CONTIGINFO | metis::METIS_DBG_CONNINFO);
        p.ptype = Some(Ptype::Pmetis);
        for ctype in [Ctype::Rm, Ctype::Shem] {
            p.ctype = Some(ctype);
            for iptype in [Iptype::Random, Iptype::Grow] {
                p.iptype = Some(iptype);
                for no2hop in [false, true] {
                    p.no2hop = no2hop;
                    for npart in [2, 4, 13, 23] {
                        p.nparts = npart;
                        p.call_suite_assert_eq();
                    }
                }
            }
        }
    }

    #[test]
    #[ignore = "slow"]
    fn little_bit_of_everthing_kmetis() {
        let mut p = Params::new(GpMetis, "", 0);
        p.nooutput = true;
        // p.dbglvl = Some(metis::METIS_DBG_IPART | metis::METIS_DBG_MOVEINFO | metis::METIS_DBG_REFINE | metis::METIS_DBG_COARSEN | metis::METIS_DBG_CONTIGINFO | metis::METIS_DBG_CONNINFO);
        p.ptype = Some(Ptype::Kmetis);
        for objtype in [Objtype::Cut, Objtype::Vol] {
            p.objtype = Some(objtype);
            for ctype in [Ctype::Rm, Ctype::Shem] {
                p.ctype = Some(ctype);
                for no2hop in [false, true] {
                    p.no2hop = no2hop;
                    for minconn in [false, true] {
                        p.minconn = minconn;
                        for contig in [false, true] {
                            p.contig = contig;
                            for npart in [2, 4, 13, 23] {
                                p.nparts = npart;
                                p.call_suite_assert_eq();
                            }
                        }
                    }
                }
            }
        }
    }
}
