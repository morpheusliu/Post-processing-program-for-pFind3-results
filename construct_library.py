file_database=open('database.txt','rt')
n=0;
while True:
	data_database=file_database.readline()
	if not data_database:
		break
	if data_database.startswith('>'):
		n+=1
		if n%1000==0:
			print(n)
		arry_data_database=data_database.split('|')
		file_id=open('library/'+arry_data_database[1],'wt')
		file_description=open('description/'+arry_data_database[1],'wt')
		file_description.write(arry_data_database[2])
	else:
		data_database=data_database.replace('\n','')
		file_id.write(data_database)
file_id.close()
file_description.close()