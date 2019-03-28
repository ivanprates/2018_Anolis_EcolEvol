GENERAL-INFO-START

		seq-file		ortonii_t70_s10_n23.gphocs
		trace-file		ortonii_n23_2000loci_400k_run1.out
		num-loci		2000
		burn-in		0
		
		mcmc-iterations		400000
		iterations-per-log		100
		logs-per-line		100

		tau-theta-print		1
		tau-theta-alpha		1
		
	  # value used and reported in the paper
		tau-theta-beta		20 
		
		# value previously reported in github 
		#tau-theta-beta		1000

		mig-rate-print		1
		mig-rate-alpha		1
		
	  # value used and reported in the paper
		mig-rate-beta		0.0000002
		
		# value previously reported in github 
		#mig-rate-beta		0.00002
		
		#start-mig		10000

		locus-mut-rate		CONST

		find-finetunes		TRUE
		find-finetunes-num-steps		100
		find-finetunes-samples-per-step		100


GENERAL-INFO-END


CURRENT-POPS-START
		
		POP-START
				name		amz
				samples		orto_LSUMZH14099 d orto_H2026 d orto_H2594 d orto_UFAC0085 d orto_BM028 d orto_BM565 d orto_LSUMZH12993 d orto_LSUMZH13904 d orto_LSUMZH14163 d orto_LSUMZH14304 d orto_MTR28130 d orto_MTR18907 d orto_MTR19241 d orto_MPEG25699 d orto_MPEG27688 d orto_MPEG29465 d orto_MPEG29325 d orto_MTR977671 d
		POP-END

		POP-START
				name		af
				samples		orto_MTR09910 d orto_MTR34436 d orto_MTR12160 d orto_MTR12291 d orto_MTR34278 d 
		POP-END

CURRENT-POPS-END


ANCESTRAL-POPS-START

		
		POP-START
				name		root
				children		amz	af
				
				# value used and reported in the paper
				tau-initial		0.003
				
				# value previously reported in github
				#tau-initial		0.002
		
		POP-END

ANCESTRAL-POPS-END


MIG-BANDS-START

		BAND-START
				source		amz
				target		af
		BAND-END

		BAND-START
				source		af
				target		amz
		BAND-END
		
MIG-BANDS-END
