GENERAL-INFO-START

		seq-file		punctatus_t70_s10_n46.gphocs
		trace-file		punctatus_n46_2000loci_400k_run1.out
		num-loci		2000
		burn-in		0
		
		mcmc-iterations		400000
		iterations-per-log		100
		logs-per-line		100

		tau-theta-print		1
		tau-theta-alpha		1
		tau-theta-beta		1000

		mig-rate-print		1
		mig-rate-alpha		1
		mig-rate-beta		0.00002
		#start-mig		10000

		locus-mut-rate		CONST

		find-finetunes		TRUE
		find-finetunes-num-steps		100
		find-finetunes-samples-per-step		100


GENERAL-INFO-END


CURRENT-POPS-START

		POP-START
				name		ceamz
				samples		punc_MTR28401 d punc_MTR28593 d punc_LSUMZH13910 d punc_MPEG29314 d punc_MPEG24758 d punc_BM288 d punc_LSUMZH14336 d punc_MPEG29943 d punc_GN71 d punc_MPEG20846 d punc_MPEG28489 d punc_MTR20798 d punc_MTR976723 d punc_MTR28048 d punc_MTR21474 d punc_MPEG22415 d punc_MTR978312 d punc_MPEG26102 d
		POP-END

		POP-START
				name		af
				samples		punc_MUFAL9635 d punc_ICST764 d punc_MTR05978 d punc_IBSPCRIB0361 d punc_MTR21545 d punc_MTR12338 d punc_JFT459 d punc_JFT773 d punc_MTR17744 d punc_MTR15267 d punc_MTRX1478 d punc_MTRX1468 d punc_MTR34227 d punc_LG1299 d punc_MTR34414 d punc_MTR12511 d
		POP-END
		
				POP-START
				name		swamz
				samples		punc_LSUMZH15476 d punc_H1907 d punc_H2546 d punc_LSUMZH14100 d punc_UNIBAN1670 d punc_H1911 d punc_MTR25584 d punc_MTR18550 d punc_MPEG21348 d punc_PJD409 d punc_LSUMZH12751 d punc_LSUMZH12577 d
		POP-END

CURRENT-POPS-END


ANCESTRAL-POPS-START

		POP-START
				name		af_ceamz
				children		ceamz		af
				tau-initial		0.001
		POP-END

		POP-START
				name		root
				children		af_ceamz		swamz
				tau-initial		0.003
		POP-END
		

ANCESTRAL-POPS-END


MIG-BANDS-START

		BAND-START
				source		ceamz
				target		af
		BAND-END

		BAND-START
				source		af
				target		ceamz
		BAND-END
		
			
		BAND-START
				source		swamz
				target		af
		BAND-END

		BAND-START
				source		af
				target		swamz
		BAND-END
		
		BAND-START
				source		swamz
				target		ceamz
		BAND-END

		BAND-START
				source		ceamz
				target		swamz
		BAND-END

MIG-BANDS-END


