include { EMMAX_VCF_GROUPPING                                      } from '../../../modules/local/gwas_vcfgroup'
include { STURCTURE_VCF_ADDID as EMMAX_VCF_ADDID                   } from '../../../modules/local/structure_vcfaddid'
include { PLINK_INDEP as EMMAX_PLINK_INDEP                         } from '../../../modules/nf-core/plink/indep/main'
include { PLINK_EXTRACT as EMMAX_PLINK_EXTRACT                     } from '../../../modules/nf-core/plink/extract/main'
include { ADMIXTURE as EMMAX_ADMIXTURE                             } from '../../../modules/nf-core/admixture/main'
include { GWAS_EMMAX_PLINK_FMAT                                    } from '../../../modules/local/gwas_emmax_plink_fmat'
include { GWAS_EMMAX_COV                                           } from '../../../modules/local/gwas_emmax_cov'
include { GWAS_EMMAX_KINSHIP                                       } from '../../../modules/local/gwas_emmax_kin'
include { GWAS_EMMAX_ASSOCIATE                                     } from '../../../modules/local/gwas_emmax_assoc'
include { GWAS_EMMAX_CMPLOT                                        } from '../../../modules/local/gwas_emmax_cmplot'
include { GWAS_EMMAX_ANNOTATION                                    } from '../../../modules/local/gwas_emmax_annotation'
include { STRUCTURE_PLOT                                           } from '../../../modules/local/structure_plot'
include { VCF2DIS                                                  } from '../../../modules/local/vcf2dis'
include { POPLDDECAY                                               } from '../../../modules/local/poplddecay'
include { GET_ENRICEMENG_GENE                                      } from '../../../modules/local/get_gwas_genelist'
include { GO_KEGG  as GO_KEGG_ANNO                                 } from '../../../modules/local/go_kegg'

workflow GWAS {
    take:
    ch_vcf_tbi_annotable
    ch_vcf_tbi_structure_vcfaddid
    ch_annoevf_emmax_anno
    main:
    
    if(params.input_group && params.groupforassociates){
        println("exist")

	ch_groupa=Channel.fromPath(params.groupforassociate).splitCsv(header:true).map{it -> tuple([group:it['group']],[groupassoc:it['groupassoc']])}
        //ch_groupa.view()
        ch_group=Channel.fromPath(params.input_group).splitCsv(header:true).map{it -> tuple([group:it['group']],[sample:it['sample']])}
        //ch_group.view()
        ch_group_info=ch_groupa.combine(ch_group,by:0)
        //ch_group_info.view()
        ch_group_split=ch_group_info.map{it-> if(it[1].groupassoc=="1"){tuple(['id':it[0].group],it[2].sample)}}.groupTuple()
        //ch_group_split.view()
        ch_all_split=ch_group_info.map{it -> tuple(['id':'popu_snp'],it[2].sample)}.groupTuple()
        //ch_all_split.view()
        ch_groups=ch_all_split.concat(ch_group_split)
        //ch_groups.view()
        ch_vcf_tbi_pca_vcfaddid=ch_vcf_tbi_annotable.groupTuple()
        //ch_vcf_tbi_pca_vcfaddid.view()
        ch_groups_vcf_emmax_groupping=ch_groups.join(ch_vcf_tbi_pca_vcfaddid)
        //ch_groups_vcf_emmax_groupping.view()
         
        ch_tree =  ch_vcf_tbi_structure_vcfaddid.map{ row -> tuple(row[0],row[2])}    
    }else{
        println("not exist")
        ch_groups_vcf_emmax_groupping=ch_vcf_tbi_structure_vcfaddid
        ch_tree = ch_vcf_tbi_structure_vcfaddid.map{ row -> tuple(row[0],row[2] )} 
    }
    VCF2DIS(ch_tree)
    ch_poplddecay = ch_tree.map{it -> tuple(it[0],params.input_group,it[1])}
    POPLDDECAY(ch_poplddecay) 
    ch_groups_vcf_emmax_groupping.view()
    //ch_groups_vcf_emmax_groupping.view()
    EMMAX_VCF_GROUPPING(
        ch_groups_vcf_emmax_groupping
    )
    // EMMAX_VCF_GROUPPING.out.vcf.view()
    EMMAX_VCF_ADDID(
        EMMAX_VCF_GROUPPING.out.vcf
    )
 
    EMMAX_PLINK_INDEP(
        EMMAX_VCF_ADDID.out.vcf
    )

    ch_vcf_tbi_prunein_emmax_plink_extract=EMMAX_VCF_ADDID.out.vcf.join(EMMAX_PLINK_INDEP.out.prunein)
    EMMAX_PLINK_EXTRACT(
        ch_vcf_tbi_prunein_emmax_plink_extract
    )

    ch_bed_bim_fam_EMMAX_k_admixture=EMMAX_PLINK_EXTRACT.out.bed.join(EMMAX_PLINK_EXTRACT.out.bim).join(EMMAX_PLINK_EXTRACT.out.fam).combine(Channel.from(2,3,4,5,6,7,8,9,10)).map{
        it -> tuple(['id':"popu_snp_${it[4]}","group":"${it[0].id}"],it[1],it[2],it[3],it[4])
    }
    //ch_bed_bim_fam_EMMAX_k_admixture.view()
    EMMAX_ADMIXTURE(
        ch_bed_bim_fam_EMMAX_k_admixture
    )

    ch_vcf_tbi_chrnum_ctgnum_emmax_fmat=EMMAX_VCF_ADDID.out.vcf.join(EMMAX_VCF_ADDID.out.ctgnum)
    //ch_vcf_tbi_chrnum_ctgnum_emmax_fmat.view()
    GWAS_EMMAX_PLINK_FMAT(
        ch_vcf_tbi_chrnum_ctgnum_emmax_fmat
    )

    ch_qfile_tpedfam_emmax_cov=EMMAX_ADMIXTURE.out.ancestry_fractions.join(EMMAX_ADMIXTURE.out.cvlog).map{
              it-> tuple(['id':"${it[0].group}"],it[1],it[2])}.groupTuple().join(
                GWAS_EMMAX_PLINK_FMAT.out.emmaxfile.map{it-> tuple(['id':"${it[0].id}"],it[2])})
    GWAS_EMMAX_COV(
        ch_qfile_tpedfam_emmax_cov
    )

    STRUCTURE_PLOT(
        GWAS_EMMAX_COV.out.cverr
    )

    ch_tped_tfam_emmax_kin=GWAS_EMMAX_PLINK_FMAT.out.emmaxfile.map{it-> tuple(['id':"${it[0].id}"],it[1],it[2],it[3],it[4])}
    GWAS_EMMAX_KINSHIP(
        ch_tped_tfam_emmax_kin
    )

    ch_emmax_assoc=GWAS_EMMAX_PLINK_FMAT.out.emmaxfile.map{it-> tuple(['id':"${it[0].id}"],it[1],it[2],it[3],it[4])}.
        join(GWAS_EMMAX_KINSHIP.out.kin).
        join(GWAS_EMMAX_COV.out.cov).combine(Channel.fromPath(params.input_phenotypes)).
        combine(Channel.fromPath(params.input_phenotypes).
        splitCsv().first().flatMap().filter{it != 'sample'}.
        merge(Channel.from(1..10000)).combine(Channel.from('K','KQ')))
   //ch_emmax_assoc.view()

    GWAS_EMMAX_ASSOCIATE(
        ch_emmax_assoc
    )
    
    GWAS_EMMAX_CMPLOT(
        GWAS_EMMAX_ASSOCIATE.out.assoc
    )

    ch_annoevf_emmax_anno=ch_annoevf_emmax_anno
    ch_signif_annoevf_emmax_anno=GWAS_EMMAX_CMPLOT.out.signif.combine(ch_annoevf_emmax_anno)
    //ch_signif_annoevf_emmax_anno.view()
    GWAS_EMMAX_ANNOTATION(
        ch_signif_annoevf_emmax_anno
    )
    if(params.species != ''){
        GET_ENRICEMENG_GENE(
        GWAS_EMMAX_ANNOTATION.out.anno_enrichment.flatten()
        )
        GO_KEGG_ANNO(GET_ENRICEMENG_GENE.out.get_gene.map{
            row -> tuple(["id":"enrichment"],row,params.species)
        })
    }
    println("GWAS workflow succeed2")
}
