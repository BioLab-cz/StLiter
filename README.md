# StLiter 81.20.0101

Release Date: 21th February, 2020

Author

	~Changyong Yu (Northeastern University in CHINA)
	~Chu Zhao (Northeastern University in CHINA)

1.Introduction
--

  StLiter is a tool，it can construct the compressed de bruijn graph (cDBG) of one genome or multiple genomes. The built graph can represent the genomes and index them in the data analysis such as large-scale read mapping.

  The k-mers of the genome sequences are the vertex of the graph, the length of the k-mer is an important parameter of the graph. StLiter limits the length of k-mer less than or equal to 128. In our experience, 128 is long enough for various kinds of aims of data analysis when using cDBG. But also, we can modify the codes to support longer k-mer vertex.

  The compressed de bruijn graph is composed of three parts of information, 1) the branching k-mers, 2) the uni-paths and 3) the postion lists of the edges. But, in our opinior, for large-scale, especially huge-scale genome or multiple genomes, we recommend FM-index for indexing the position of the edges for reducing the cost of memory. 

  StLiter is a flexible tool. You can use StLiter to design the best workflow for constructing cDBG with some k-mer length according to the features and sizes of the genomes. The other advantage of cDBG is that it can process genomes of very large-scale, in our experiments it can deal with up to 100 human-size genomes with a pc of 160 GB memory on acceptable time cost.

  There is a shortcoming of StLiter. StLiter is designed base on the assumption that the alphbet is {A,C,G,T}. Therefore, StLiter can just recognize the four letters. It cannot deal with the other characters with recognizing them. If it encounter with other characters, it will treat them as letter 'A'. This surely will result in false positive results. The further verification on the genome sequence will correct and delete the false positive results. It is worth attention for the users of cDBG for this shortcoming. We think it is not vitally important compared with the key problem of how to filter out true negitive results in data analysis such as read mapping. Therefore, it is accecptable.

  Overall, StLiter can effectively construct cDBG for multiple large-scale genomes for read mapping. We believe that the tool will improve the performance of data analysis with the bottleneck of large-scale data.

2.Test Data
--

  StLiter constructs compressed de bruijn graph of genome sequences. The genome sequences are .fa format, and we tested StLiter with the genome sequences used in paper "TwoPaCo: An efficient algorithm to build the compressed de Bruijn graph from many complete genomes". we download the data using website: https://github.com/medvedevgroup/TwoPaCo.

3.Building Notes
--

Our directory contains only a content : 1.codes to build compressed de bruijn graph.

The directory named #src# consisted of the code necessary to build the StLiter.
To build the StLiter, change the directory to src and type

	$ make

After that, you can type commands to run tests using StLiter. StLiter provides 5 methods for building  compressed de bruijn graph of large-scale genomes.

4.Usage Notes
--

1)B method. Generate branching k-mers of the cDBG with short k-mer length.

	$ ./StLiter -M <method> -L <kmer_len> -t <thread_num> -r <ref_path> -c <rc_label> -a <adjacent_label> -s [single_task_length]
	
  B method calcultes  effectively the branching k-mers when 'L' is small (L<=log(4,Maximal memory)). It output two files, take L=16 for example, it will output two files, "kmer016in" and "kmer016out" respectively. They are the files saving the in-branching and out-branching k-mers with k=16. The k-mers in the files are sorted with the ascending lexicographic order. 

Running test:

	$ ./StLiter -M B -L 16 -t 4 -r dataset1 -c 1 -a 1

  The above commond takes "dataset1" file and "nomLeu2.fa" file as input, it will output kmer016in, kmer016inad, kmer016out, kmer016outad kmer016uni and kmer016uniad binary files. Each branching k-mer is saved as an uint64_t type data in the file. Adjacent information of each k-mer is saved as an uint8_t type data. The commond costs at least 4^16 Bytes memory.

2)E method. Extend short branching k-mers to long branching k-mers.

	$ ./StLiter -M <method> -L <kmer_len> -t <thread_num> -r <ref_path> -l <extend_len> -i <in_or_out> -a <adjacent_label> -p [position_vector_label] -s [single_task_length] -S [sort_result_label] -m <hashtable_memory> -c <rc_label> 

  E method calculates branching k-mers of long k-mer length based on the branching k-mers of short k-mer length. 

Running test:

	$ ./StLiter -M E -L 16 -t 4 -r dataset1 -l 3 -i 0 -a 1 -p 3 -m 10 -c 1
		
	Input files:
		kmer016in, dataset1, nomLeu2.fa and posi016in000

	Output files:
		kmer019in,kmer019inad and posi019in000

  The above commond calculates in-branching k-mers of length 19 based on the given branching k-mers with length 16. The extend length is 3.

3)F method. Extend short branching k-mers to long branching k-mers.

	$ ./StLiter -M <method> -L <kmer_len> -t <thread_num> -r <ref_path> -l <extend_len> -i <in_or_out> -a <adjacent_label> -p [position_vector_label] -s [single_task_length] -S [sort_result_label] -m <hashtable_memory> -c <rc_label> 

  F method calculates branching k-mers of long k-mer length based on the branching k-mers of short k-mer length. Different with E method, F method deals with the sparse cases of candidates. For example, if 'extend_len=10', then for each original branching k-mer, 4^10 candidates are checked. However, it is possible that only a few of the candidates, such as ~1000 candidates, are present actually. In this case, E method may cost much memory meaninglessly. F method is more effective for solving problem in this case. Generally, when the length of the k-mer is bigger than some value, then it is turn to using F method instead of E method. For "nomLeu2.fa" data, when the length of the k-mer is bigger than 35, we use F method other than E method. 

Running test:

	$ ./StLiter -M F -L 35 -t 4 -r dataset1 -l 10 -i 0 -a 1 -p 3 -m 10 -c 1
		
	Input files:
		kmer035in, dataset1, nomLeu2.fa and posi035in000

	Output files:
		kmer045in,kmer045inad and posi045in000

  The above commonds calculate branching k-mers of length 45 based on the given branching k-mers with length 35. The extend length is 10.

4)U method. Merge in-branching and out-branching k-mers to branching k-mers.

	 $ ./StLiter -M <method> -L <kmer_len> -a <adjacent_label>

  U method merges in-branching and out-branching k-mers to branching k-mers.

Running test:

	$ ./StLiter -M U -L 35 -a 1

	Input files:
		kmer035out, kmer035in, kmer035inad and kmer035outad

	Output files:
		kmer035,kmer035ad

  The above commonds merge kmer035out and kmer035in into kmer035, and merger kmer035inad and kmer035outad into kmer035ad. 

5)P method. Calculate the uni-paths and position list of edges.

	$ ./StLiter -M <method> -L <kmer_len> -t <thread_num> -r <ref_path> -p [position_vector_label] -c <rc_label> -m <hashtable_memory>

  P method calcultes uni-paths and save each unipath as an uint64_t type data. The first 39 bits saves the starting position of the uni-path, the next 17 bits saves the length of the uni-path and the last 8 bits saves the genome id. 

Running test:

	$ ./StLiter -M P -L 35 -t 4 -r dataset1 -p 3 -c 1 -m 5 -P 0

	Input files:
		kmer035, dataset1, nomLeu2.fa, posi035in000 and posi035out000

	Output files:
		unipath035 and counterfile

  The above commonds calculate the unipaths and position lists.

5.Parameter Settings
--

The format of a parameter of StLiter in the command line is a pair of strings, here we denote the pair as (-p, [q]) or (-p,<q>). String p is the name of the parameter. String q is the value of the parameter input in the command line. [q] represents that the parameter is a optional parameter. <q> represents that the parameter is a necessary parameter.

@parameter (-r,<ref_path>)

  Parameter 'r' gives the path of a text file which saves the filenames of the genome files. For example 'ref_path'="dataset1", than dataset1 is a text file which saves the filename "nomLeu2.fa".

@parameter (-L,<kmer_len>)

  Note that parameter 'L' not only set the length of the k-mer, but also decide the memory cost. The program will take at least 4^kmer_len Bytes memory (lm). That is,

	k=15	lm=1GB
	k=16	lm=4GB
	k=17	lm=16GB
	k=18	lm=64GB

@parameter (-s,[single_task_length])

  The parameter 's' is optional, the default value for that parameter is 24, which means that, in the multiple thread computation each calculating thread will process 2^24 length of sequence at each time of task. When dealing with short genomes, the parameter should be sett smaller than the number of genome_len/thread_num. 

@parameter (-i,<in_or_out>)

  Parameter 'i' assigns the case of the computation, 'in_or_out=0' indicates the case of in-branching k-mer calculations, 'in_or_out=1' indicates the case of out-branching k-mer calculation.

@parameter (-a,<adjacent_label>)

  Parameter 'a' decides whether the adjacent information of the branching k-mers is saved and output. 'adjacent_label=0' indicates the case of un-saving adjacent information of branching k-mer. 'adjacent_label=1' indicates saving adjacent information.

@parameter (-l,<extend_len>)

  Parameter 'l' assigns the extending length in this step, for example if 'extend_len=5', the original k-mer length is 'L', then the k-mer length will become 'L+5' after extending process. It is noting that parameter 'l' actually decides the total number of possible branching k-mers checked in the current step. The total number is N*4^l, where N is the total number of branching k-mer with length 'L'. According to our experience, we recommend that parameter 'l' should be set satisfying that N*4^l<64G.

@parameter (-p,[position_vector_label])

  Parameter 'p' decides whether using the bit-vector which labels the possible starting positions of branching k-mers along the genome sequence. The length of bit-vector is equal to that of labeled genome sequence. 'position_vector_label=0' denotes un-using bit-vector and un-saving bit-vector of current step. 'position_vector_label=1' denotes un-using bit-vector and saving bit-vector of current step. 'position_vector_label=2' denotes using bit-vector and un-saving bit-vector of current step. 'position_vector_label=3' denotes using bit-vector and saving bit-vector of current step.
	Noting that in P method, parameter 'p' just labels whether using bit-vector. 'position_vector_label=2' or 'position_vector_label=3' indicate using bit-vector, otherwise un-use. P method can not generate new bit-vectors.

@parameter (-S,[sort_result_label])

  Parameter 'S' decides whether sorting the branching k-mer results. 'sort_result_label=1' means sorting the results and 'sort_result_label=0' means un-sorting the results.

@parameter (-m, <hashtable_memory>)

  Parameter 'm' decides the maximal memory cost of the hash table used by E or F method in GB. Note that the memory cost of the hash table is not the total memory cost. Actually, the total memory cost of StLiter is hard to be estimated precisely with different genome sequences. In our experience, the parameter can be set as follows.

	Memory available	hashtable_memory	
	>=32GB			5
	>=80GB			10
	>=160GB			15
	other			20

@parameter (-c, <rc_label>)

  Parameter 'c' decide whether calculate the reverse complement strand of the genome sequence. 'rc_label=1' indicates calculating the CDBG of the genome sequences and their reverse complement sequences. 'rc_label=0' indicates calculating the CDBG of only the genome sequences without their reverse complements.

@parameter (-P, <position_list_label>)

  Parameter 'P' decide whether calculating and output the position list of the edges. 'position_list_label=1' indicates calculating and output the position list of the edges. 'rc_label=0' indicates not calculating and output the position list of the edges.


6.Format of cDBG
--

	-kmerxxxin
	-kmerxxxout
	-kmerxxx
		~binary files of branching k-mers
		~xxx is the k-mer length
		~each r Bytes as a k-mer,where r=8*([kmer_len/32]+sign(kmer_len/32)).

	-kmerxxxinad
	-kmerxxxoutad
	-kmerxxxad
		~binary files of adjacent information of branching k-mers
		~xxx is the k-mer length
		~the i-th Byte assigns the adjacent information of the i-th k-mer in the corresponding files. 

	-posixxxinyyy
	-posixxxoutyyy
		~binary file of bit vector of genomes
		~xxx is the k-mer length
		~yyy is the id of the genome
		~the i-th bit label the branching information of the i-th position in the genome with id:yyy.
		
	-unipathxxx
		~binary file of uni-path
		~xxx is the k-mer length
		~each 8 Bytes as a uni-path, the first 39 bits indicate the starting position, the next 17 bits indicate the length and the last 8 bits indicate the genome id. 

	-poslistxxx
		~binary file of position lists
		~xxx is the k-mer length
		~the position list of each edge is saved in three parts 
		1)first r Bytes as a (k+1)-mer, 
		2)the next 8 Bytes as the length of the list and 
		3) the positions, each position of which is saved in 8 Bytes. where r=8*([(kmer_len+1)/32]+sign((kmer_len+1)/32)). The position lists of all edges are saved one by one. 	

7.License
--

	See LICENSE.txt

8.Contacts
--

	Please e-mail your feedback at cyyneu@126.com.



