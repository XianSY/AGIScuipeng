include { TWAS_fun } from '../../../modules/local/twas'

workflow TWAS{
    
    take:
        twas_param

    main:
        //test_param = [["id":"group"],[ch_expression],[ch_GM],[ch_pheno]]
        //println(test_param)
	TWAS_fun(
            twas_param
        )

}
