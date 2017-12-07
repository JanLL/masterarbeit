import numpy as np
import matplotlib.pyplot as plt



if __name__ == '__main__':


	filename = 'output.txt'

	T_list = []
	cp_list = []
	dcp_list = []

	with open(filename) as f:
		for line in f:

			#line = line.replace('(a)', '')
			#line = line.replace('\n', '')
			line = line.split('\t')
			line = [float(x) for x in line]

			T_list.append(line[0])
			cp_list.append(line[1])
			dcp_list.append(line[2])


	fig = plt.figure()
	ax = fig.add_subplot(111)
	#ax.plot(T_list, cp_list, label='c_p')
	#ax.plot(T_list, dcp_list, label='dc_p')
	ax.plot(T_list, np.array(cp_list) - np.array(dcp_list))
	ax.set_xlabel('T [degC]')
	ax.legend(loc='upper left')
	
	
	plt.show()

#	for row in data:




