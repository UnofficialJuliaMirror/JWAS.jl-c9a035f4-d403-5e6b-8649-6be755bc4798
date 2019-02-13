################################################################################
# Pre-Check
################################################################################
function errors_args(mme,methods)
    Pi         = mme.MCMCinfo.Pi
    estimatePi = mme.MCMCinfo.estimatePi
    if methods == "conventional (no markers)"
        if mme.M!=0
            error("Conventional analysis runs without genotypes!")
        elseif estimatePi == true
            error("conventional (no markers) analysis runs with estimatePi = false.")
        end
    elseif methods=="RR-BLUP"
        if mme.M == 0
            error("RR-BLUP runs with genotypes")
        elseif Pi != 0.0
            error("RR-BLUP runs with π=0.")
        elseif estimatePi == true
            error("RR-BLUP runs with estimatePi = false.")
        end
    elseif methods=="BayesC"
        if mme.M == 0
            error("BayesC runs with genotypes.")
        end
    elseif methods=="BayesB"
        if mme.M==0
            error("BayesB runs with genotypes.")
        end
    elseif methods=="BayesL"
        if mme.M == 0
            error("BayesL runs with genotypes.")
        elseif estimatePi == true
            error("BayesL runs with estimatePi = false.")
        end
    elseif methods=="GBLUP"
        if mme.M == 0
            error("GBLUP runs with genotypes.")
        elseif mme.M.genetic_variance == false
            error("Please provide values for the genetic variance for GBLUP analysis")
        end
    end
    if mme.nModels > 1 && mme.M!=0
        if Pi != 0.0 && round(sum(values(Pi)),digits=2)!=1.0
          error("Summation of probabilities of Pi is not equal to one.")
        end
    end
end

function check_pedigree(mme,df,pedigree)
    if mme.ped == 0 && pedigree == false
        return
    end
    if pedigree!=false
        pedID=map(string,collect(keys(pedigree.idMap)))
    else
        pedID=map(string,collect(keys(mme.ped.idMap)))
    end

    if mme.M!=0 && !issubset(mme.M.obsID,pedID)
        error("Not all genotyped individuals are found in pedigree!")
    end

    phenoID = map(string,df[1])
    if !issubset(phenoID,pedID)
        error("Not all phenotyped individuals are found in pedigree!")
    end
end

function check_outputID(mme)
    #Genotyped individuals are usaully not many, and are used in GWAS (complete
    #and incomplete), thus are used as default output_ID if not provided
    if mme.M == 0 && mme.pedTrmVec == 0
        mme.MCMCinfo.outputEBV = false #no EBV in non-genetic analysis
    end

    if mme.MCMCinfo.outputEBV == false
        mme.output_ID = 0
    elseif mme.output_ID == 0 && mme.M != 0
        mme.output_ID = mme.M.obsID
    elseif mme.output_ID == 0 && mme.M == 0 && mme.pedTrmVec != 0
        #output EBV for all individuals in the pedigree for PBLUP
        pedID=map(string,collect(keys(mme.ped.idMap)))
        mme.output_ID = pedID
    end
end

function check_phenotypes(mme,df)
    single_step_analysis = mme.MCMCinfo.single_step_analysis
    phenoID = map(string,df[1])   #same to df[:,1] in deprecated CSV
    if mme.M == 0 && mme.ped == 0 #non-genetic analysis
        return df
    end
    if single_step_analysis == false && mme.M != 0 #complete genomic data
        if !issubset(phenoID,mme.M.obsID)
            printstyled("Phenotyped individuals are not a subset of\n",
            "genotyped individuals (complete genomic data,non-single-step).\n",
            "Only use phenotype information for genotyped individuals.\n",bold=false,color=:red)
            index = [phenoID[i] in mme.M.obsID for i=1:length(phenoID)]
            df    = df[index,:]
            printstyled("The number of individuals with both genotypes and phenotypes used\n",
            "in the analysis is ",size(df,1),".\n",bold=false,color=:red)
        elseif mme.output_ID!=0 && !issubset(mme.output_ID,mme.M.obsID)
            printstyled("Testing individuals are not a subset of \n",
            "genotyped individuals (complete genomic data,non-single-step).\n",
            "Only output EBV for tesing individuals with genotypes.\n",bold=false,color=:red)
            mme.output_ID = intersect(mme.output_ID,mme.M.obsID)
        end
    else #incomplete genomic data , PBLUP
        pedID = map(string,collect(keys(mme.ped.idMap)))
        if !issubset(phenoID,pedID)
            printstyled("Phenotyped individuals are not a subset of\n",
            "individuals in pedigree (incomplete genomic data (single-step) or PBLUP).\n",
            "Only use phenotype information for individuals in the pedigree.\n",bold=false,color=:red)
            index = [phenoID[i] in pedID for i=1:length(phenoID)]
            df    = df[index,:]
            printstyled("The number of individuals with both phenotype and pedigree information\n",
            "used in the analysis is ",size(df,1),".\n",bold=false,color=:red)
        elseif mme.output_ID!=0 && !issubset(mme.output_ID,pedID)
            printstyled("Testing individuals are not a subset of \n",
            "individuals in pedigree (incomplete genomic data (single-step) or PBLUP).\n",
            "Only output EBV for tesing individuals in the pedigree.\n",bold=false,color=:red)
            mme.output_ID = intersect(mme.output_ID,pedID)
        end
    end
    return df
end

function init_mixed_model_equations(mme,df,sol)
    getMME(mme,df)
    #starting value for sol can be provided
    if sol == false #no starting values
        sol = zeros(T,size(mme.mmeLhs,1))
    else            #besure type is Float64
        sol = map(T,sol)
    end
    return sol,df
end

function speedup(mme,starting_value)
    starting_value = map(Float32,starting_value)
    mme.X          = map(Float32,mme.X)
    mme.ySparse    = map(Float32,mme.ySparse)
    mme.mmeLhs     = map(Float32,mme.mmeLhs)
    mme.mmeRhs     = map(Float32,mme.mmeRhs)
    if mme.pedTrmVec != 0 #if no SSBR
        mme.Ai = map(Float32,mme.Ai)
    end
    return starting_value
end
