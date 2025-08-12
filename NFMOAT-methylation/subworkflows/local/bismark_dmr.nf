include { BISMARK2BEDGRAPH as BISMARK2BEDGRAPH_CHG    } from '../../modules/local/bismark2graph'
include { BISMARK2BEDGRAPH as BISMARK2BEDGRAPH_CHG    } from '../../modules/local/bismark2graph'
include { BISMARK2BEDGRAPH as BISMARK2BEDGRAPH_CHH    } from '../../modules/local/bismark2graph'
include { METILENE_INPUT   as METILENE_INPUT_CPG      } from '../../modules/local/metilene_input'
include { METILENE_INPUT   as METILENE_INPUT_CHG      } from '../../modules/local/metilene_input'
include { METILENE_INPUT   as METILENE_INPUT_CHH      } from '../../modules/local/metilene_input'
include { METILENEDMR                                 } from '../../modules/local/metilenedmr'

workflow BISMARK_DMR{
    take:
    methylation_chg
    methylation_cpg
    methylation_chh

    main:
    ch_bismark2bedgraph_CHG = BISMARK_METHYLATIONEXTRACTOR.out.methylation_chg.map{
        it -> tuple(it[0],it[1],it[2])
    }
    ch_bismark2bedgraph_CpG = BISMARK_METHYLATIONEXTRACTOR.out.methylation_cpg.map{
        it -> tuple(it[0],it[1],it[2])
    }
    ch_bismark2bedgraph_CHH = BISMARK_METHYLATIONEXTRACTOR.out.methylation_chh.map{
        it -> tuple(it[0],it[1],it[2])
    }
    /*
     * 分开获取得到的CpG、CHH、CHG的 beggraph 文件
     * 1、首先将CpG、CHH、CHG的calling 文件分开成三个输入文件
    */
    ch_group = channel.fromPath(params.group).map{it -> tuple(['id':'group'],it)}
    BISMARK2BEDGRAPH_CHG(
        ch_bismark2bedgraph_CHG
    )   
    BISMARK2BEDGRAPH_CPG(
        ch_bismark2bedgraph_CpG
    )
    BISMARK2BEDGRAPH_CHH(
        ch_bismark2bedgraph_CHH
    )
    //BISMARK2BEDGRAPH_CHG.out.bedgraph.map{it -> it[1]}.collect().view()
    
    METILENE_INPUT_CHG(   
    ch_group,
    BISMARK2BEDGRAPH_CHG.out.bedgraph.map{it -> it[1]}.collect()
    )
    METILENE_INPUT_CPG(
    ch_group,
    BISMARK2BEDGRAPH_CPG.out.bedgraph.map{it -> it[1]}.collect()
    )  
    METILENE_INPUT_CHH(
    ch_group,
    BISMARK2BEDGRAPH_CHH.out.bedgraph.map{it -> it[1]}.collect()
    )
    ch_compare = channel.fromPath(params.compare).splitCsv(header:true).map{row -> tuple(row.formerGroup,row.latergroup)}
    ch_dmrbedgraph = METILENE_INPUT_CHG.out.group_bedgraph.join(METILENE_INPUT_CpG.out.group_bedgraph).join(METILENE_INPUT_CHH.out.group_bedgraph).groupTuple()
    ch_dmrbedgraph.view()
    METILENEDMR(
     ch_compare,
     ch_dmrbedgraph      
    )
    output:

}