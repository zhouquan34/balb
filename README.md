# Copyright 1/16/2016 Quan Zhou, Baylor College of Medicine 
# Joint work with Philip Ernst and L. C. G. Rogers 

'balb' is for simulation of Bermudan/American Look-Back option using method in [1]. 
License: GPL-3

You can simply use the binary executables 'balb-linux' or 'balb-mac' depending on your operating system. 

To view the Help, type:

	./balb-mac help

An example command:
	
	./balb-mac -a 5 -o test -t 200

This will simulate the American look-back option expiring on day 200 with a 5-day look-back window. For more details, please refer to [2]. 
The results will be saved to test.high.out and test.low.out, which contain the stopping day and reward of each simulated path.

References:
1. L. C. G. Rogers. Bermudan options by simulation. Technical report, University of Cambridge, 2015.
2. P. Ernst, L. C. G. Rogers, Q. Zhou. The value of foresight. Submitted in 2016.

