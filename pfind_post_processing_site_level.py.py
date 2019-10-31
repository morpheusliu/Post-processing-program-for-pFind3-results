calculation=input('Please enter the calculation style (lh for light/heavy & hl for heavy/light):')
intscore=input('Please enter the cut off value of interference score:')
mc=input("Please enter the number of modification in each spectrum(a for single modification & b for multi-modification):")
print(calculation)
import glob
file_list=glob.glob('data/*')
file_counts=0
for file_name in file_list:
	file=open(file_name,'rt')	#这一个模块是进行score_interference的基本筛选
	modification=file_name[file_name.find('\\')+1:file_name.find('_')]
	print(modification)
	data_quality=open('int/1_data1.txt','wt')
	file_counts+=1
	print(file_counts)
	while True:
		raw_data=file.readline()
		if not raw_data:
			break
		arry_raw_data=raw_data.split('\t')
		try:
			arry_raw_data[10]=round(float(arry_raw_data[10]),5)	#尝试将第11列的数据（score_interference）转换为浮点数字格式，这是为了排除第一行的标题行
		except:
			continue
		if arry_raw_data[10]<=float(intscore):	#interference score的筛选值 根据需要进行修改
			data_quality.write(arry_raw_data[1]+'\t'+arry_raw_data[2]+'\t'+arry_raw_data[6]+'\t'+str(round(float(arry_raw_data[9]),5))+'\t'+str(arry_raw_data[10])+'\t'+str(round(float(arry_raw_data[4]),5))+'\n')
	data_quality.close()
	file.close()
	
	file_data_quality=open('int/1_data1.txt','rt')	#这个模块是进行数据的最终筛选，去掉con_项，去掉rev_
	file_final_data=open('int/2_data2.txt','wt')
	while True:
		data1=file_data_quality.readline()
		if not data1:
			break
		arry_data1=data1.split('\t')
		arry2_arry_data1=arry_data1[2].split('$')
		xcount=0
		for acs in arry2_arry_data1:
			if 'CON_' in acs:	#计数con_项，但是可能有些项不是污染但是含有con_，在后两行将这种情况筛除
				xcount=xcount+1
			if '|' in acs and 'CON_' in acs:	#因为污染项不含有|，所以用这个排除可能含有con但是非污染的项，综合这一行和上两行的计数结果，如果xcount还是大于0则说明确实存在污染，则将这一行数据去掉
				xcount=xcount-1
		if mc=="a":
			if arry_data1[1].count(modification)==1:	#选择有且仅有一个修饰关键词的项。这一行可以按照要求进行修改。
				if xcount<=0:
					tem=arry_data1[1][:arry_data1[1].find(','+modification)]
					site=tem[tem.rfind('|')+1:]
					arry_arry_data1=arry_data1[2].split('$')
					arry_all_ac=[]
					for ac in arry_arry_data1:
						if '|' in ac:
							ac1=ac[ac.find('|')+1:ac.rfind('|')]
							arry_all_ac.append(ac1)
					if len(arry_all_ac)>0:	#筛选掉去掉反向库匹配后没有匹配到任何一个有效ac的二级质谱数据项
						print(arry_data1[0]+'\t'+site+'\t'+(' ').join(arry_all_ac)+'\t'+arry_data1[3]+'\t'+arry_data1[4]+'\t'+arry_data1[5],file=file_final_data,end='')
		if mc=="b":
			arry_sites=[]
			if modification in arry_data1[1]:
				if xcount<=0:
					modification_inf=arry_data1[1]
					while modification_inf:
						if not modification in modification_inf:
							break
						modification_inf=modification_inf[:modification_inf.rfind(','+modification)]
						arry_sites.append(modification_inf[modification_inf.rfind('|')+1:])
					arry_arry_data1=arry_data1[2].split('$')
					arry_all_ac=[]
					for ac in arry_arry_data1:
						if '|' in ac:
							ac1=ac[ac.find('|')+1:ac.rfind('|')]
							arry_all_ac.append(ac1)
					if len(arry_all_ac)>0:	#筛选掉去掉反向库匹配后没有匹配到任何一个有效ac的二级质谱数据项
						print(arry_data1[0]+'\t'+(' ').join(arry_sites)+'\t'+(' ').join(arry_all_ac)+'\t'+arry_data1[3]+'\t'+arry_data1[4]+'\t'+arry_data1[5],file=file_final_data,end='')
	file_data_quality.close()
	file_final_data.close()
	
	print('final data complete')
	print(str(file_counts)+'Protein data is filtering')
	
	import os
	import shutil
	shutil.rmtree('data_protein_group')
	os.mkdir('data_protein_group')
	
	file_protein_data=open('data_protein/'+file_name[file_name.find('\\')+1:]+'_p','rt')
	file_names=''
	types='a'
	while True:
		data_protein=file_protein_data.readline().replace('\n','')
		arry_data_protein=data_protein.split('\t')
		if not data_protein:
			break
		if len(arry_data_protein)>=2:
			try:
				int(arry_data_protein[0])
			except:
				if types=='b':
					if 'SameSet'==arry_data_protein[1]:
						if '|' in arry_data_protein[2]:
							nn+=1
							if nn<=20:	#将seameset的数量限制在20个以内
								file_result=open('data_protein_group/'+file_names,'at')
								print(arry_data_protein[2][arry_data_protein[2].find('|')+1:arry_data_protein[2].rfind('|')],file=file_result,end=' ')
								file_result.close()
				continue
			if int(arry_data_protein[0])>0 and '|' in arry_data_protein[1]:
				file_names=arry_data_protein[1][arry_data_protein[1].find('|')+1:arry_data_protein[1].rfind('|')]
				file_result=open('data_protein_group/'+file_names,'wt')
				file_result.close()
				nn=0
				types='b'
	file_protein_data.close()
	
	file_final_data=open('int/2_data2.txt','rt')	#这一个模块用于匹配修饰肽段的位点在对应蛋白ac上的位点位置（这是个核心模块，一般不进行修改）。19年6月12日进行了修改，原来位点匹配只进行一次，现在leading和sameset分别进行位点匹配
	file_match_site=open('int/3_confirm_site.txt','wt')
	n=0
	while True:
		final_data=file_final_data.readline()
		if not final_data:
			break
		arry_final_data=final_data.split('\t')
		arry_final_data[0]=arry_final_data[0].replace('I','L')
		arry_tem2=arry_final_data[2].split(' ')
		arry_sitess=arry_final_data[1].split(' ')
		arry_ac_site_leading=[]
		nnn=0
		for acsssss in arry_tem2:
			if os.path.exists('data_protein_group/'+acsssss):
				nnn+=1
		if nnn>1:	#进行unique peptide筛选
			continue
		arry_all_site=[]
		for ac in arry_tem2:	#第一次匹配，这次进行leading的位点匹配
			if os.path.exists('data_protein_group/'+ac):
				file_library=open('library/'+ac,'rt')
				while True:
					sequence=file_library.readline()
					break
				file_library.close()
				sequence=sequence.replace('I','L')
				if sequence.count(arry_final_data[0])>1:
					continue
				for sitess in arry_sitess:
					while True:
						site=sequence.rfind(arry_final_data[0])
						if not site+1:
							break
						final_site=site+int(sitess)
						arry_all_site.append(ac+'_'+str(final_site)+'_'+str(sitess))
						sequence=sequence[:site]
		for acb in arry_all_site:	#leading匹配完了之后，根据leading匹配的结果，找到leading在pfind鉴定中的所有sameset，再次进行位点匹配
			acb_final=acb[:acb.find('_')]
			file_acs=open('data_protein_group/'+acb_final)	#这里确定leading的sameset是从pfind的group结果里直接把sameset结果调进来，sameset已经在group收集的时候进行了筛选，全都是正常的ac
			while True:
				data_acs=file_acs.readline().replace('\n','')
				break
			file_acs.close()
			arry_data_acs=data_acs.split(' ')
			arry_sameset_ac=[]
			for acss in arry_data_acs:
				if acss:
					file_library=open('library/'+acss,'rt')
					while True:
						sequence=file_library.readline()
						break
					file_library.close()
					sequence=sequence.replace('I','L')
					if sequence.count(arry_final_data[0])>1:
						continue
					while True:
						site=sequence.rfind(arry_final_data[0])
						if not site+1:
							break
						final_site=site+int(acb[acb.rfind('_')+1:])
						sequence=sequence[:site]
						arry_sameset_ac.append(acss+'_'+str(final_site))
			print(final_data.replace('\n','')+'\t'+acb[:acb.rfind('_')]+'\t'+(' ').join(arry_sameset_ac),file=file_match_site)
		n+=1
		if n%500==0:
			print('adding sites:processing data'+str(n))
	file_final_data.close()
	file_match_site.close()
	
	print('Site confirmation complete')
	print(str(file_counts)+' constructing leader ac list')
	
	file_data3=open('int/3_confirm_site.txt','rt')	#这个模块是把文件3里的数据完整的拷贝到文件4中，而文件4用于后续的分析工作。
	file_data4=open('int/4_format_match_site.txt','wt')
	while True:
		data_data3=file_data3.readline()
		if not data_data3:
			break
		print(data_data3,file=file_data4,end='')
	file_data3.close()
	file_data4.close()
	
	file_list=open('int/5_list.txt','wt')
	file_list.close()
	
	file_format_match_site=open('int/3_confirm_site.txt','rt')	#这个模块用于从site confirmation结果找出所有涉及的leading ac 并收集其对应的sequence和site
	arry_leader_ac=[]
	while True:
		data_format_match_site=file_format_match_site.readline().replace('\n','')
		if not data_format_match_site:
			break
		arry_data_format_match_site=data_format_match_site.split('\t')
		arry_arry6_format_match_site=arry_data_format_match_site[6].split(' ')
		for ac_site in arry_arry6_format_match_site:
			if not ac_site in arry_leader_ac and ac_site:
					arry_leader_ac.append(ac_site)
	file_format_match_site.close()
	for leader_ac in arry_leader_ac:
		arry_all_seq=[]
		file_format_match_site=open('int/4_format_match_site.txt','rt')
		arry_file_format_match_site=file_format_match_site.readlines()
		file_format_match_site.close()
		if not arry_file_format_match_site:	#当文档读完则跳出循环（因为这个文档在运行过程中内容不断减少，最后减少到零）
			break
		file_data5=open('int/4_format_match_site.txt','wt')
		for ac_x in arry_file_format_match_site:
			data_format_match_site=ac_x.replace('\n','')
			if not data_format_match_site:
				break
			if leader_ac in data_format_match_site:
				arry_data_format_match_site=data_format_match_site.split('\t')
				arry_arry6_format_match_site=arry_data_format_match_site[6].split(' ')
				if leader_ac in arry_arry6_format_match_site:
					arry4_data_format_match_site=arry_data_format_match_site[1].split(' ')
					final_seq=arry_data_format_match_site[0]
					for sitesss in arry4_data_format_match_site:
						final_seq=final_seq[:int(sitesss)]+'('+modification+')'+final_seq[int(sitesss):]	#把原始序列修改为添加了修饰的序列
					if not final_seq in arry_all_seq:
						arry_all_seq.append(final_seq)
			else:
				print(ac_x,file=file_data5,end='')
		file_data5.close()
		file_list=open('int/5_list.txt','at')
		print(leader_ac[:leader_ac.find('_')]+'\t'+leader_ac[leader_ac.find('_')+1:]+'\t'+(' ').join(arry_all_seq),file=file_list)
		file_list.close()
	
	file_data3=open('int/3_confirm_site.txt','rt')	#这个模块是把文件3里的数据完整的拷贝到文件4中，而文件4用于后续的分析工作。
	file_data15=open('int/15_temp1.txt','wt')
	while True:
		data_data3=file_data3.readline()
		if not data_data3:
			break
		print(data_data3,file=file_data15,end='')
	file_data3.close()
	file_data15.close()
	
	file_list_2=open('int/5_list.txt','rt')	#这个模块根据list的leading结果重新回到site confirm结果 找到所有对应的sameset结果，由于sameset还附带了位点信息，无法直接从group结果中找sameset
	file_list_3=open('int/6_final_list.txt','wt')
	while True:
		data_list_2=file_list_2.readline().replace('\n','')
		if not data_list_2:
			break
		arry_list_2=data_list_2.split('\t')
		arry_sameset=[]
		file_data11=open('int/15_temp1.txt','rt')
		arry_data11=file_data11.readlines()
		file_data11.close()
		if not arry_data11:
			break
		file_data12=open('int/15_temp1.txt','wt')
		for data11 in arry_data11:
			data_final=data11.replace('\n','')
			if '\t'+arry_list_2[0]+'_'+arry_list_2[1] in data_final:
				arry_data_final=data_final.split('\t')
				arry_arry6_data_final=arry_data_final[6].split(' ')
				if arry_list_2[0]+'_'+arry_list_2[1] in arry_arry6_data_final:
					arry_arry7_data_final=arry_data_final[7].split(' ')
					for ac_sites in arry_arry7_data_final:
						if not ac_sites in arry_sameset:
							arry_sameset.append(ac_sites)
			else:
				print(data11,file=file_data12,end='')
		file_data12.close()
		print(data_list_2+'\t'+(' ').join(arry_sameset)+'\t',file=file_list_3)
	file_list_2.close()
	file_list_3.close()

	print('final list is complete')
	print(str(file_counts)+' calculating the test result')
	
	file_data3=open('int/3_confirm_site.txt','rt')	#这个模块是把文件3里的数据完整的拷贝到文件4中，而文件4用于后续的分析工作。
	file_data16=open('int/16_temp2.txt','wt')
	while True:
		data_data3=file_data3.readline()
		if not data_data3:
			break
		print(data_data3,file=file_data16,end='')
	file_data3.close()
	file_data16.close()
	
	file_final_list=open('int/6_final_list.txt','rt')
	file_calculating_result=open('int/7_calculating_result.txt','wt')
	while True:
		data_final_list=file_final_list.readline()
		arry_final_list=data_final_list.split('\t')
		arry_ratio=[]
		arry_interference=[]
		arry_identification=[]
		arry_file_match_site_plus=[]
		arry_all_sequ=[]
		file_match_site=open('int/16_temp2.txt','rt')
		arry_file_match_site_plus=file_match_site.readlines()
		file_match_site.close()
		if not arry_file_match_site_plus:
			break
		file_datax=open('int/16_temp2.txt','wt')
		for ac_xx in arry_file_match_site_plus:
			data_match_site=ac_xx.replace('\n','')
			arry_data_match_site=data_match_site.split('\t')
			if arry_final_list[0]+'_'+arry_final_list[1]==arry_data_match_site[6]:
				arry_ratio.append(arry_data_match_site[3])
				arry_interference.append(arry_data_match_site[4])
				arry_identification.append(arry_data_match_site[5])
			else:
				print(ac_xx,file=file_datax,end='')
		file_datax.close()
		if calculation=='hl':
			sum_ratio=0
			for rat in arry_ratio:
				sum_ratio+=float(rat)
			ave_ratio=sum_ratio/len(arry_ratio)
			sum_ratio_sd=0
			for rat in arry_ratio:
				sum_ratio_sd+=(float(rat)-ave_ratio)**2
			if len(arry_ratio)>1:
				ratio_sd=(sum_ratio_sd/(len(arry_ratio)-1))**0.5
			elif len(arry_ratio)==1:
				ratio_sd='NA'
		if calculation=='lh':
			sum_ratio=0
			for rat in arry_ratio:
				sum_ratio+=1/float(rat)
			ave_ratio=sum_ratio/len(arry_ratio)
			sum_ratio_sd=0
			for rat in arry_ratio:
				sum_ratio_sd+=(1/float(rat)-ave_ratio)**2
			if len(arry_ratio)>1:
				ratio_sd=(sum_ratio_sd/(len(arry_ratio)-1))**0.5
			elif len(arry_ratio)==1:
				ratio_sd='NA'
		if ratio_sd=='NA':
			ratio_sd='NA'
		else:
			ratio_sd=round(ratio_sd,5)
		sum_inte=0
		for inte in arry_interference:
			sum_inte+=float(inte)
		ave_inte=sum_inte/len(arry_ratio)
		sum_inte_sd=0
		for inte in arry_interference:
			sum_inte_sd+=(float(inte)-ave_inte)**2
		if len(arry_ratio)>1:
			inte_sd=(sum_inte_sd/(len(arry_ratio)-1))**0.5
		elif len(arry_ratio)==1:
			inte_sd='NA'
		sum_ide=0
		for ide in arry_identification:
			sum_ide+=float(ide)
		ave_ide=sum_ide/len(arry_ratio)
		sum_ide_sd=0
		for ide in arry_identification:
			sum_ide_sd+=(float(ide)-ave_ide)**2
		if len(arry_ratio)>1:
			ide_sd=(sum_ide_sd/(len(arry_ratio)-1))**0.5
		elif len(arry_ratio)==1:
			ide_sd='NA'
		if inte_sd=='NA':
			inte_sd='NA'
		else:
			inte_sd=round(inte_sd,5)
		if ide_sd=='NA':
			ide_sd='NA'
		else:
			ide_sd=round(ide_sd,5)
		print(data_final_list.replace('\n','')+str(round(ave_ratio,5))+'\t'+str(ratio_sd)+'\t'+str(round(ave_inte,5))+'\t'+str(inte_sd)+'\t'+str(round(ave_ide,5))+'\t'+str(ide_sd)+'\t'+str(len(arry_ratio)),file=file_calculating_result)
	file_final_list.close()
	file_calculating_result.close()
	
	print('Calculation is complete')
	print(str(file_counts)+' Adding description')
	
	file_calculating_result=open('int/7_calculating_result.txt','rt')
	file_description=open('int/8_add_description.txt','wt')
	while True:
		data_calculating_result=file_calculating_result.readline()
		if not data_calculating_result:
			break
		arry_data_calculating_result=data_calculating_result.split('\t')
		file_description_lib=open('description/'+arry_data_calculating_result[0],'rt')
		while True:
			data_description_lib=file_description_lib.readline()
			break
		file_description_lib.close()
		print(data_calculating_result.replace('\n','')+'\t'+data_description_lib,file=file_description,end='')
	file_calculating_result.close()
	file_description.close()
	
	print('Description complete')
	
	print('Description complete')
	print(str(file_counts)+' Adding genename')
	
	file_description=open('int/8_add_description.txt','rt')
	file_genename=open('int/9_add_genename.txt','wt')
	while True:
		data_description=file_description.readline().replace('\n','')
		if not data_description:
			break
		arry_data_description=data_description.split('\t')
		if 'GN=' in arry_data_description[11]:
			genename_l=arry_data_description[11][arry_data_description[11].find('GN='):]
			genename_l=genename_l[:genename_l.find(' ')]
			genename_l=genename_l.replace('GN=','')
		else:
			genename_l=''
		arry_sames=arry_data_description[3].split(' ')
		arry_all_genename=[]
		for ss in arry_sames:
			try:
				file_description_lib=open('description/'+ss[:ss.find('_')],'rt')
			except:
				continue
			while True:
				data_description_lib=file_description_lib.readline()
				break
			file_description_lib.close()
			gnc=0
			for gns in arry_all_genename:
				if 'GN='+gns+' ' in data_description_lib:
					gnc+=1
			if gnc==0:
				try:
					data_description_lib=data_description_lib[data_description_lib.find('GN='):]
					data_description_lib=data_description_lib[:data_description_lib.find(' ')]
					data_description_lib=data_description_lib.replace('GN=','')
				except:
					continue
				if data_description_lib:
					arry_all_genename.append(data_description_lib)
		print(data_description+'\t'+genename_l+'\t'+(' ').join(arry_all_genename),file=file_genename)
	file_description.close()
	file_genename.close()
	
	file_genename=open('int/9_add_genename.txt','rt')
	file_result=open('results/site_'+calculation+'_'+file_name.replace('data\\',''),'wt')
	print('AC\tSite\tGenename\tSameSets\tAll Genename_SameSets\tModified sequences\tSpectum Count\tDescription\tRatio '+calculation+'\tSD Ratio\tInterference Score\tSD int\tIdentification Score\tSD ide',file=file_result)
	while True:
		data_genename=file_genename.readline().replace('\n','')
		if not data_genename:
			break
		arry_data_genename=data_genename.split('\t')
		print(str(arry_data_genename[0])+'\t'+str(arry_data_genename[1])+'\t'+str(arry_data_genename[12])+'\t'+str(arry_data_genename[3])+'\t'+str(arry_data_genename[13])+'\t'+str(arry_data_genename[2])+'\t'+str(arry_data_genename[10])+'\t'+str(arry_data_genename[11])+'\t'+str(arry_data_genename[4])+'\t'+str(arry_data_genename[5])+'\t'+str(arry_data_genename[6])+'\t'+str(arry_data_genename[7])+'\t'+str(arry_data_genename[8])+'\t'+str(arry_data_genename[9]),file=file_result)
	file_genename.close()
	file_result.close()
	
	file_data10=open('int/9_add_genename.txt','rt')	#这个模块用于收集leading和sameset的所有ac用于功能注释分析，作为数据的下载列表
	file_list_ac=open('results/list_'+calculation+'_'+file_name.replace('data\\',''),'wt')
	while True:
		data_data10=file_data10.readline().replace('\n','')
		if not data_data10:
			break
		arry_data10=data_data10.split('\t')
		print(arry_data10[0],file=file_list_ac)
		arry_arry3_data10=arry_data10[3].split(' ')
		for ac_s in arry_arry3_data10:
			acsb=ac_s[:ac_s.find('_')]
			print(acsb,file=file_list_ac)
	file_data10.close()
	file_list_ac.close()
	
	print(str(file_counts)+' complete')