# Responsive threshold search based memetic algorithm for balanced minimum sum-of-squares clustering (BMSSC)
This repository includes the source code of the proposed MA algorithm published in an Information Sciences paper titled with "Responsive threshold search based memetic algorithm for balanced minimum sum-of-squares clustering".

The 16 datasets used are available from the UCI machine learning repository (http://www.ics.uci.edu/mlearn/MLRepository.html). To facilitate the further research, we upload the instances here.

We made comparisons between MA and some state-of-the-art methods from the following related BMSSC works.

[1] L.R. Costa, D. Aloise, N. Mladenovic´ , Less is more: basic variable neighborhood search heuristic for balanced minimum sum-of-squares clustering, Information Sciences 415–416 (2017) 247–253.

[2] M.I. Malinen, P. Fränti, Balanced K-means for clustering, in: P. Fränti, G. Brown, M. Loog, F. Escolano, M. Pelillo (Eds.), Structural, Syntactic and Statistical Pattern Recognition, S+SSPR 2014. Lecture Notes in Computer Science, vol. 8621, Springer, Berlin, 2014, pp. 32–41..

Please cite our work as:

Zhou, Q., Hao, J. K., & Wu, Q. (2021). Responsive threshold search based memetic algorithm for balanced minimum sum-of-squares clustering. Information Sciences, 569, 184-204.

** Instructions to use the source code of MA

*** To compile:

q.zhou$ make

q.zhou$

*** To run:

q.zhou$ ./MA_BMSSC ./input_file n d k time seed ./output_res_file

(where input_file is the instance name, n is the numbe of points of the instance, d is the dimension of the point, k is the number of clusterrs, 
		time is the cutoff of an execution, seed is the random seed, such as 1, 2, ..., 10., output_res_file is a file used to store the running information)

q.zhou$

*** To clean

q.zhou$ make clean

q.zhou$
