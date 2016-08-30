#!/usr/bin/python
#############################################################################
# Genome Assembly and Annotation Error Estimation Algorithm Utilizing CAFEv3.0
# Written by: Gregg Thomas, Summer 2012
# Hahn Lab, Indiana University
# Contact: grthomas@indiana.edu
#
# This version of the script first minimizes scores for varying error models across all species 
# and then (if specified) incrementally adds or subtracts 10% of the global error to each species 
# to minimize scores and estimate error for each species individually.
#
# May 2016: Added an error message that displays if a lambda command with more than one search
# parameter is input. We don't recommend doing this. We do recommend estimating error with a
# 1 lambdda model and then subsequently using that error on a multi-lambda run of CAFE.
# Streamlined the output code with the printWrite function.
# Removed getAllSpec function and replaced it with regex.
#############################################################################

import sys, argparse, os, random, datetime, time, re

############################################
#Function Definitions
############################################
def errorOut(errnum, errmsg):
# Formatting for error messages.
	fullmsg = "|**Error " + str(errnum) + ": " + errmsg + " |";
	border = " " + "-" * (len(fullmsg)-2);
	print "\n" + border + "\n" + fullmsg + "\n" + border + "\n";

############################################
def optParse(errorflag):
# This function handles the command line options.

	parser = argparse.ArgumentParser();

	parser.add_argument("-i", dest="input_file", help="A CAFE shell script with the full CAFE path in the shebang line, the load, tree, and lambda commands. These lines will be read and incorporated into the caferror shell script.");
	parser.add_argument("-e", dest="user_err_start", help="The starting point for the grid search. Should be between 0 and 1. Default: 0.4", type=float, default=0.4);
	parser.add_argument("-d", dest="user_tmp_dir", help="A directory in which all caferror files will be stored. If none is specified, it will default to caferror_X, with X being some integer one higher than the last directory.", default="");
	parser.add_argument("-f", dest="first_run", help = "Boolean option to perform a pre-error model run (1) or not (0). Default: 0", type=int, default=1);
	parser.add_argument("-c", dest="curve_option", help="Boolean option. caferror can either perform the grid search (0) or search a pre-specified space (1). Default: 0", type=int, default=0);
	parser.add_argument("-t", dest="error_tries", help="A list of error values to search over. Note: -c MUST be set to 1 to use these values. Enter as a comma delimited string, ie -t 0.1,0.2,0.3", default="")
	parser.add_argument("-l", dest="user_log_file", help="Specify the name for caferror's log file here. Default: caferrorLog.txt", default="caferrorLog.txt");
	parser.add_argument("-o", dest="output_file", help="Output file which stores only the error model and score for each run. Default: caferror_default_output.txt", default="caferror_default_output.txt");
	parser.add_argument("-s", dest="ind_min", help="Boolean option to specify whether to perform only the global error search (0) or continue with individual species minimizations (1). Default: 0", type=int, default=0);
	parser.add_argument("-v", dest="verbose", help="Boolean option to have detailed information for each CAFE run printed to the screen (1) or not (0). Default: 1", type=int, default=1);
	parser.add_argument("-m", dest="run_mode", help=argparse.SUPPRESS, type=int, default=0);

	args = parser.parse_args();

	if errorflag == 0:

		if args.input_file == None:
			errorOut(1, "-i must be defined");
			optParse(1);

		if not all(op in [0,1] for op  in [args.first_run, args.curve_option, args.ind_min, args.verbose]):
			errorOut(2, "-f, -c, -s, and -v must all take values of either 0 or 1");
			optParse(1);

		if args.user_err_start > 1 or args.user_err_start < 0:
			errorOut(3, "-e must take values of between 0 and 1");
			optParse(1);

		if args.curve_option == 0 and args.error_tries != "":
			errorOut(4, "With -t specified -c must also be set to 1");
			optParse(1);

		if args.user_tmp_dir != "" and args.user_tmp_dir[len(args.user_tmp_dir) - 1] != "/":
			args.user_tmp_dir = args.user_tmp_dir + "/";

		if args.curve_option == 1 and args.error_tries == "":
			args.error_tries = [0.0, 0.001, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95];
		elif args.curve_option == 1 and args.error_tries != "":
			args.error_tries = args.error_tries.split(",");
		
		### -m: Run mode for individual species minimzations (for debugging purposes only... not supported).
		###		0 - Default: No shuffle, background constant.
		###		1 - No shuffle, background updated.
		###		2 - Shuffle, background constant.
		###		3 - Shuffle, background updated.

		return args.input_file, args.output_file, args.user_err_start, args.user_tmp_dir, args.curve_option, args.error_tries, args.user_log_file, args.ind_min, args.run_mode, args.verbose, args.first_run;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

############################################
def printWrite(o_name, v, o_line1, o_line2="", pad=0):
# Function to print a string AND write it to the file.
	if o_line2 == "":
		outline = o_line1;
	else:
		outline = o_line1 + " "*(pad-len(o_line1)) + o_line2;
	if v == 1 or v == -2:
		print outline;
	f = open(o_name, "a");
	f.write(outline + "\n");
	f.close();

############################################
def cafeRun(families, newick, lamb, path, tempDir, cafLog, specDict, spectomin, logrm, fflag, nr, vb):
#This is the function which creates an error file based on the current error (by calling maxfamsize), writes the appropriate caferror.sh script
#and executes the current caferror.sh script.

	printWrite(cafLog, 1, "--------------------------");

	outline = "CAFE run number " + str(nr);
	if spectomin[0] == "final":
		outline = outline + "(FINAL RUN)";
	printWrite(cafLog, 1, outline);

	if spectomin[0] == "all":
		printWrite(cafLog, 1, str(specDict[errSpec[0]]) + " Error Model Run for all");
		logFile = tempDir + "cafe_" + str(specDict[errSpec[0]]) + "_all_log.txt";
	elif spectomin[0] == "none":
		printWrite(cafLog, 1, "Pre run (no error models)");
		logFile = tempDir + "cafe_pre_log.txt";
	elif spectomin[0] == "final":
		printWrite(cafLog, 1, "With minimized error models");
		logFile = tempDir + "cafe_final_log.txt";
		repFile = tempDir + "cafe_final_report";
	else:
		printWrite(cafLog, 1, str(specDict[spectomin[0]]) + " Error Model Run for " +  spectomin[0]);
		logFile = tempDir + "cafe_" + str(specDict[spectomin[0]]) + "_" + spectomin[0] + "_log.txt";
	printWrite(cafLog, 0, "-----");
	# The above logic statements name the log file for the current CAFE run. The format for the log files is: cafe_[errormodel]_[species]_log.txt.
	# Special cases for species include 'all' in which the current errormodel is applied to all during global error estimation species and 'final' 
	# which is the last CAFE run with the minimized errormodels applied to each individual species.

	if logrm == 1:
		os.system("rm " + logFile);
	#In the event that CAFE initializes incorrectly, it must be re-run with the same parameters. When that happens, the old log file is removed here
	#and will be replaced with one of the same name. The contents of the logfile are saved to the 'badRuns.txt' file in the getScore function.

	if len(spectomin) == 1:
	#For global error prediction (spectomin = ['all']), the final CAFE run (spectomin = ['final']), and individual species minimizations (spectomin = [the
	#current species] the spectomin list will contain just one element.
		if spectomin[0] == "none":
			printWrite(cafLog, 0, "## No error models generated");	
		else:
			for key in specDict:
				ErrString = "cafe_errormodel_" + str(specDict[key]) + ".txt";
				ErrFile = tempDir + ErrString

				genErrFile = errFileCheck(key, ErrString, specDict, tempDir, cafLog)
				#This line calls the errFileCheck function to see if the error model file for this error distribution and this species is already created.

				if genErrFile == 1:
					maxfamsize(families, specDict[key], ErrFile);
					#This line calls the maxfamsize function below to create the error file.

	# else:
	# #When -g is specified in the command line, spectomin will contain more than one element. [NOTE: -g not yet developed]
	# 	for key in spectomin:

	# 		ErrString = "cafe_errormodel_" + str(specDict[key]) + ".txt";
	# 		ErrFile = tempDir + ErrString
		
	# 		genErrFile = errFileCheck(key, ErrString, specDict, tempDir, cafLog)
	# 		#This line calls the errFileCheck function to see if the error model file for this error distribution and this species is already created.

	# 		if genErrFile == 1:
	# 			maxfamsize(families, specDict[key], ErrFile);
	# 			#This line calls the maxfamsize function below to create the error file.
	# # The above blocks create the error file for the current run (if necessary) by calling the errFileCheck and maxfamsize functions.

	#####
	printWrite(cafLog, 0, "-----");
	printWrite(cafLog, 1, "Rewriting CAFE shell script...");

	##########
	# These lines write the caferror.sh script.
	shellfile = "caferror.sh";

	# shellfileList = os.listdir(path);
	# shellcount = 1;
	# while shellfile in shellfileList:
	#	shellfile = "caferror" + str(shellcount) + ".sh";
	#	shellcount = shellcount + 1;
	# cLog.write("\nshellfile: " + shellfile + "\n");
	# cLog.write("logfile: " + logFile + "\n");
	# Some lines for debugging. Uncomment as desired.

	cafeFile = open(shellfile, "w");

	cafeFile.write(path + "\n" + newick + "\n");

	if fflag == 0:
		loadLine = "load -i " + families + " -t 10 -l " + logFile + "\n"
	elif fflag == 1:
		loadLine = "load -i " + families + " -t 10 -l " + logFile + " -filter\n"

	cafeFile.write(loadLine);

	if spectomin[0] != "none":
		if len(spectomin) == 1:
			for key in specDict:
				ErrString = "cafe_errormodel_" + str(specDict[key]) + ".txt";
				ErrFile = tempDir + ErrString;
				errLine = "errormodel -model " + ErrFile + " -sp " + key + "\n";
				cafeFile.write(errLine);	
	
		else:
			for key in spectomin:
				ErrString = "cafe_errormodel_" + str(specDict[key]) + ".txt";
				ErrFile = tempDir + ErrString;
				errLine = "errormodel -model " + ErrFile + " -sp " + key + "\n";
				cafeFile.write(errLine);	

	cafeFile.write(lamb + "\n");
	if spectomin[0] == "final":
		repline = "report " + repFile;
		cafeFile.write(repline);

	cafeFile.close();
	# End caferror.sh writing.
	##########

	os.system("chmod +x " + shellfile);
	if vb == 0:
		printWrite(cafLog, 1, "Running CAFE [silently] with error models listed above...");
		os.system("./" + shellfile + " >> " + tempDir + "cafe.out 2>&1");
	else:
		printWrite(cafLog, 1, "Running CAFE with error models listed above...");
		os.system("./" + shellfile);
	nr = nr + 1;
	#These lines give the caferror.sh script permissions and execute it to run cafe given the current parameters.

	print("CAFE run complete! Retrieving Score..........");
	#This is printed here. Score retrieval is actually handled by getScore (directly below), which is always called immediately after 
	#this function in the main block.

	return nr;

############################################
def errFileCheck(errspec_check, errstring_check, specdict_check, thedir, cafLog):
#This function checks if, for a given species and error model, the error model file has already been created. It returns 1 the error model file
#has not yet been created and 0 if it already has.

	fileList = os.listdir(thedir);

	if errstring_check in fileList:
		printWrite(cafLog, 0, "## " + str(specdict_check[errspec_check]) + " Error File Already Created for " + errspec_check);
		return 0;

	else:
		printWrite(cafLog, 0, "## Generating " + str(specdict_check[errspec_check]) + " Error File for " + errspec_check);	
		return 1;

############################################
def getScore(error, spectomin, tempDir, wout, cafLog):
#This function retrieves the score calculated by CAFE from a given log file, specified by the amount of error and species used in that run.
#This function now also checks if the previous CAFE run was initialized properly. If so, the program will continue and if not it will send
#a signal back to the call telling it to re-run CAFE with the same parameters.

	init = "good";
	lcount = 0;
	linfcount = 0;

	if spectomin[0] == "all":
		scoreFile = tempDir + "cafe_" + str(error) + "_all_log.txt";
	elif spectomin[0] == "final":
		scoreFile = tempDir + "cafe_final_log.txt";
	elif spectomin[0] == "none":
		scoreFile = tempDir + "cafe_pre_log.txt";
	else:
		scoreFile = tempDir + "cafe_" + str(error) + "_" + spectomin[0] + "_log.txt";

	sFile = open(scoreFile, "r");
	sLines = sFile.readlines();
	sFile.close();

	initFile = open(tempDir + "InitFile.txt", "a");
	initFile.write(scoreFile);
	initFile.write("\n")

	for s in xrange(len(sLines)):
		if sLines[s][:7] == "Poisson":
			initFile.write(sLines[s]);
		#if sLines[s][:7] == "Poisson" and sLines[s].find("inf") != -1:
		#	init = "bad";
		if sLines[s][:7] == ".Lambda":
			lcount = lcount + 1;
		if sLines[s][:7] == ".Lambda" and sLines[s].find("inf") != -1:
			linfcount = linfcount + 1;

		if sLines[s].find("Lambda Search Result:") != -1:
			lamval = sLines[s+1][9:sLines[s+1].index("&") - 1];
			score = sLines[s+1][sLines[s+1].index("Score: ") + 7:].replace("\n","");
			#score = score.replace('\n', '');

	if lcount == linfcount:
		init = "bad";

	initFile.write(str(score));
	initFile.write("\n");
	initFile.write(init);
	initFile.write("\n\n");
	initFile.close();


	printWrite(cafLog, 1, "Score with above error models:", str(score), 35);
	printWrite(cafLog, 1, "Lambda with above error models:", str(lamval), 35);


	if init == "bad":
		printWrite(cafLog, 1, "++WARNING: CAFE failed to converge or initialize. This run will be re-done. Check badRuns.txt for more info.\n");

		bFile = open(tempDir + "badRuns.txt", "a");
		bFile.write("\n************\n\n")
		bFile.write(scoreFile);
		bFile.write("\n");
		for sline in sLines:
			bFile.write(sline);
		bFile.close();

	if init == "good":
		if wout == 1:
			outFile = open(outFilename, "a");
			line = str(error) + "\t" + score + "\n"
			outFile.write(line);
			outFile.close();

	return float(score), lamval, init;

############################################
def maxfamsize(inFilename, totError, errFilename):
# This function creates an error file in the proper format for CAFE by reading the input gene family file and
# finding the gene family with the largest number of genes.

	negAsym = 0.5;
	posAsym = 1 - negAsym;
	# I used this to run some simulations with asymmetric error distributions. Right now it is set to split the error evenly
	# between +1 and -1, but go ahead and change negAsym if needed.

	inFile = open(inFilename, "r");
	lines = inFile.readlines();
	inFile.close();

	maxfs = 0;

	for i in range(len(lines)):
		if i == 0:
			continue;		
		cline = lines[i].replace("\n","").split("\t");	
		k = 2;
	
		while k <= len(cline) - 1:
			if int(cline[k]) > maxfs:
				maxfs = int(cline[k]);			
			k = k + 1;

	erroutFile = open(errFilename, "w");
	erroutFile.write("maxcnt:");
	erroutFile.write(str(maxfs));
	erroutFile.write("\n");
	erroutFile.write("cntdiff -1 0 1\n");

	j = 0;

	while j <= maxfs:
		if j == 0:
			pline = str(j) + " 0.00 " + str((totError / 2) + (1 - totError)) + " " + str(totError / 2);	
			erroutFile.write(pline);
		else:
			pline = str(j) + " " + str(totError * negAsym) + " " + str(1 - totError) + " " + str(totError * posAsym);		
			erroutFile.write(pline);

		if j != maxfs:
			erroutFile.write("\n");
		j = j + 1;
	
	erroutFile.close();

############################################
#Main Block
############################################

startsec = time.time();
start = datetime.datetime.now().strftime("%m.%d.%Y | %I:%M:%S");

inFilename, outFilename, initError, tmpDir, wholeCurveOpt, errTries, caferrorLog, indSpecMin, Mode, vOpt, firstRun = optParse(0);
# The first step is to retrieve the proper values for options and filenames based on the user's command line specifications.
# errTries = [0.0, 0.001, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95];

initCheck = "";
spec_to_min = [];
rmlog = 0;
FilterFlag = 0;

inFile = open(inFilename, "r");
inLines = inFile.readlines();
inFile.close();

for x in range(len(inLines)):
	if inLines[x][:2] == "#!":
		CafePath = inLines[x].replace("\n","");

	if inLines[x][:4] == "load":
		FamFile = inLines[x][(inLines[x].index("-i") + 3):];
		FamFile = FamFile[:FamFile.index(" ")];

		if inLines[x].find("-filter") != -1:
			FilterFlag = 1;

	if inLines[x][:4] == "tree":
		Tree = inLines[x].replace('\n', '');

	if inLines[x][:6] == "lambda":
		LamStruct = inLines[x].replace('\n', '');

		if "-t" in LamStruct:
			t_struct = LamStruct[LamStruct.index("-t"):];
			if not all(c in '-t \t(1,);' for c in t_struct):
				errorOut(5, "Only single lambda searches are recommended for estimating error. For more info, see the CAFE Manual.");
				sys.exit();
# Obtains input variables from an existing CAFE script specified with -i.

treestring = Tree[Tree.index("("):];
treestring = re.sub('[)][\d.eE-]+:[\d.eE-]+', ')', treestring);
treestring = re.sub(':[\d.eE-]+', '', treestring);
treestring = re.sub('<[\d]+>', '', treestring);
errSpec = treestring.replace("(","").replace(")","").replace(";","").split(",");
# This block extracts the species names from the input phylogeny.

mainSpecDict = {};
for eachspec in errSpec:
	mainSpecDict[eachspec] = initError;
spec_to_min.append("all");
# Some initializations of important variables above. 
# errSpec gets a list of all species in the input tree. 
# Throughout the program, mainSpecDict will keep track of the species and their corresponding minimized errors. It follows the [key]:[value] format
# of [species]:[minimized error model].
# spec_to_min is what will be passed as the species to minimize for the global error estimation algorithm. Usually this will be passed simply 
# as ['all'].

#########
# This block creates a new default directory for each run of caferror.py
if tmpDir == "":
	prevList = os.listdir(os.getcwd());

	isDirMade = 0;
	x = 1

	while isDirMade == 0:
		dirName = "caferror_" + str(x);
		if dirName not in prevList:		
			tmpDir = os.getcwd() + "/" + dirName + 	"/";
			isDirMade = 1;
		else:
			x = x + 1;

	if tmpDir.find("//") != -1:
		tmpDir = tmpDir.replace("//", "/");

	os.system("mkdir " + tmpDir);
	# If the user didn't define a temp directory with -d, the default value is set here. This also handles a possible error with the double slashes.

else:
	os.system("mkdir " + tmpDir);

# End default dir block
##########

caferrorLog = tmpDir + caferrorLog;
if caferrorLog.find("//") != -1:
	caferrorLog = caferrorLog.replace("//", "/");

clogfile = open(caferrorLog, "w");
clogfile.write("");
clogfile.close();
# If the user didn't define a log file with -l, the default value is set here. This also handles a possible error with the double slashes.

outFilename = tmpDir + outFilename;
outFile = open(outFilename, "w");
outFile.write("ErrorModel\tScore\n");
outFile.close();
# This primes the output file with the column headers.

####################################
### Begin input info block!
pad = 35

printWrite(caferrorLog, 1, "# =========================================================================");
printWrite(caferrorLog, 1, "#\t\tAssembly/annotation error estimation");
printWrite(caferrorLog, 1, "#\t\t\t" + str(start));
printWrite(caferrorLog, 1, "# Using CAFE shell file:", inFilename, pad);
printWrite(caferrorLog, 1, "# --------------------------");
printWrite(caferrorLog, 1, "#\t\t\tINPUT INFO");
printWrite(caferrorLog, 1, "# CAFE path set as:", CafePath, pad);
printWrite(caferrorLog, 1, "# Using gene family file:", FamFile, pad);
printWrite(caferrorLog, 1, "# Using tree command:", Tree, pad);
printWrite(caferrorLog, 1, "# Using lambda command", LamStruct, pad);
printWrite(caferrorLog, 1, "# --------------------------");
printWrite(caferrorLog, 1, "#\t\t\tOPTIONS INFO");
printWrite(caferrorLog, 1, "# -c " + str(wholeCurveOpt), "Global grid search option.", pad);
printWrite(caferrorLog, 1, "# -e " + str(initError), "Starting error estimate of " + str(initError) + ".", pad);
printWrite(caferrorLog, 1, "# -f " + str(firstRun), "Initial CAFE run option.", pad);
printWrite(caferrorLog, 1, "# -v " + str(vOpt), "CAFE verbosity option.", pad);
printWrite(caferrorLog, 1, "# -s " + str(indSpecMin), "Individual species minimization option.", pad);
printWrite(caferrorLog, 1, "# --------------------------");
printWrite(caferrorLog, 1, "#\t\t\tOUTPUT INFO");
printWrite(caferrorLog, 1, "# Putting all files in:", tmpDir, pad);
printWrite(caferrorLog, 1, "# Log file for all runs:", caferrorLog, pad);
printWrite(caferrorLog, 1, "# Output file:", outFilename, pad);
printWrite(caferrorLog, 1, "# =========================================================================");
### End input info block!
####################################

printWrite(caferrorLog, 1, "# Beginning Global Error Prediction...");

numruns = 1;

if firstRun == 1:
	spec_to_min[0] = "none";
	while initCheck == "bad" or initCheck == "":
		numruns = cafeRun(FamFile, Tree, LamStruct, CafePath, tmpDir, caferrorLog, mainSpecDict, spec_to_min, rmlog, FilterFlag, numruns, vOpt);
		lastScore, lastLam, initCheck = getScore(initError, spec_to_min, tmpDir, 1, caferrorLog);
		if initCheck == "bad":
			rmlog = 1;
	initCheck = "";
	rmlog = 0;
	preScore = lastScore;
	preLam = lastLam;
	spec_to_min[0] = "all";
#The pre-error model run if -f 1.

####################################
#Global error estimation begins here for -c 0 (default)

if wholeCurveOpt == 0:
#If -c is set to 0, which is the default value, caferror enters this block to predict error in an iterative fashion that uses information from the 
#previous run until a minimum score is achieved. Caferror starts at the error value specified by -e (0.4 by default).
	errList = [];
	errList.append(initError);

	errMin = initError;
	#Initializations of some variables. The 'Holder' variables are used to keep the place of parameters which are un-needed for any one cafeRun call.

	while initCheck == "bad" or initCheck == "":
		numruns = cafeRun(FamFile, Tree, LamStruct, CafePath, tmpDir, caferrorLog, mainSpecDict, spec_to_min, rmlog, FilterFlag, numruns, vOpt);
		lastScore, lastLam, initCheck = getScore(initError, spec_to_min, tmpDir, 1, caferrorLog);
		if initCheck == "bad":
			rmlog = 1;
	initCheck = "";
	rmlog = 0;
	#The initial Cafe run using the error valued specified by -e (0.4 by default).

	minScore = lastScore;
	errMin = initError;
	lastError = initError;

	nextError = 0.0;
	posLimit = 1.0;
	posLimScore = lastScore;
	negLimit = 0.0;
	negLimScore = 0.0; 

	currentError = initError / 2;
	errList.append(currentError);

	for eachspec in mainSpecDict:
		mainSpecDict[eachspec] = currentError;

	tally = 0;
	keepGoing = 1;

	while keepGoing == 1:
	#This will keep guessing error models until one of two termination scenarios are met. These scenarios are checked for at the
	#end of each CAFE run below.

		tally = tally + 1;

		errList.append(currentError);
		while initCheck == "bad" or initCheck == "":
			numruns = cafeRun(FamFile, Tree, LamStruct, CafePath, tmpDir, caferrorLog, mainSpecDict, spec_to_min, rmlog, FilterFlag, numruns, vOpt);
			currentScore, currentLam, initCheck = getScore(currentError, spec_to_min, tmpDir, 1, caferrorLog);
			if initCheck == "bad":
				rmlog = 1;
		initCheck = "";
		rmlog = 0;

#		outFile = open(outFilename, "a");
#		line = str(currentError) + "\t" + str(currentScore) + "\n"
#		outFile.write(line);
#		outFile.close();
#		print("++++++++++";
#		print errMin, minScore;
#		print negLimit, posLimit;
#		print("++++++++++";
#		cLog.write("errMin and minScore before logic: ");
#		cLog.write(str(errMin));
#		cLog.write(", ");
#		cLog.write(str(minScore));
#		cLog.write("\n");
#		cLog.write("negLimit and posLimit before logic: ");
#		cLog.write(str(negLimit));
#		cLog.write(", ");
#		cLog.write(str(posLimit));
#		cLog.write("\n");
		#Some print statements to track values for debugging. Uncomment as desired.

		##################
		#The main logic statements of the code are below. These take the score and error values of the most recent CAFE run
		#and evaluate them against some information obtained from previous runs, eventually arriving at a point 
		#where the termination scenarios are met with a minimum score value.
		if currentScore < minScore:
			if currentError < lastError:
				nextError = currentError / 2;

				posLimit = lastError * 2;
				posLimScore = lastScore;

			elif currentError > lastError:
				nextError = currentError * 2;

				negLimit = lastError / 2;
				megLimScore = lastScore;

			minScore = currentScore;	
			minLam = currentLam;
			errMin = currentError;

		elif currentScore > minScore:
			if currentError < errMin:
				nextError = errMin * 2

				negLimit = currentError;
				negLimScore = currentScore;

			elif currentError > errMin:
				nextError = errMin / 2;

				posLimit = currentError;
				posLimScore = currentScore;

#		print("----------";
#		print errMin, minScore;
#		print negLimit, posLimit;
#		print("----------";
#		cLog.write("errMin and minScore after logic: ");
#		cLog.write(str(errMin));
#		cLog.write(", ");
#		cLog.write(str(minScore));
#		cLog.write("\n");
#		cLog.write("negLimit and posLimit after logic: ");
#		cLog.write(str(negLimit));
#		cLog.write(", ");
#		cLog.write(str(posLimit));
#		cLog.write("\n");
		#Some print statements to track values for debugging. Uncomment as desired.
	
		if nextError > posLimit:
			nextError = ((currentError + posLimit) / 2);
		if nextError < negLimit:
			nextError = ((currentError + negLimit) / 2);
		if nextError in errList:
			nextError = ((currentError + nextError) / 2);

#		print("***************";
#		print currentError;
#		print len(str(currentError));
#		print("***************";
		#Some print statements to track values for debugging. Uncomment as desired.

		lastError = currentError;

		#####
		#These are the termination scenarios. If the error parameter reaches a certain level of precision (12 decimal places), or if
		#CAFE has been run 20 times, error estimation is terminated. This may be improved.
		if len(str(currentError)) >= 15:
			keepGoing = 0;
		elif tally > 20:
			keepGoing = 0;
		#End termination block.
		#####
		else:
			currentError = nextError;
			for eachspec in mainSpecDict:
				mainSpecDict[eachspec] = nextError;
			lastScore = currentScore;
		#End logic block.
		##################

#Global error estimation ends here for -c 0 (default)
####################################



####################################
#Global error estimation begins here for -c 1

elif wholeCurveOpt == 1:
#If -c is set to 1 by the user caferror enters this block to calculate scores for a pre-set set of values [0,1).
#This will be less precise than the method above, but will give the user an easy to visualize output when plotted.

	minScore = 1000000000.0;
	errMin = 1.0;
	#Initializations of some variables.


	for currentError in errTries:
	#This runs Cafe using all error models in the list 'errTries' defined above. Feel free to add or remove values as desired.

		for eachspec in mainSpecDict:
			mainSpecDict[eachspec] = currentError;
	
		while initCheck == "bad" or initCheck == "":
			numruns = cafeRun(FamFile, Tree, LamStruct, CafePath, tmpDir, caferrorLog, mainSpecDict, spec_to_min, rmlog, FilterFlag, numruns, vOpt);
			currentScore, currentLam, initCheck = getScore(currentError, spec_to_min, tmpDir, 1, caferrorLog);
			if initCheck == "bad":
				rmlog = 1;
		initCheck = "";
		rmlog = 0;

		if float(currentScore) < minScore and minScore != float("infinity"):
			minScore = float(currentScore);
			minLam = currentLam;
			errMin = currentError;


#Global error estimation ends here for -c 1
####################################


####################################
#Individual species error estimation begins here if -s 1
#
#After the above algorithm has minimized error scores across all species (by applying the same error model to all species simultaneously), the user has
#the option to continue minimization for each individual species. This is done by iteratively adding and subtracting 10% of the global min error from each
#current species's error model until the score ceases to improve.

#errMin = 0.04765625;
#minScore = 50384.197476;
#minLam = 0.27626649636157;

finalSpecDict = {};

if Mode == 2 or Mode == 3:
	random.shuffle(errSpec);

if Mode == 1 or Mode == 3:
	for eachspec in mainSpecDict:
		mainSpecDict[eachspec] = errMin;
	overallMinScore = minScore;

for eachspec in mainSpecDict:
	finalSpecDict[eachspec] = 0.0

if indSpecMin == 1:
	printWrite(caferrorLog, 1, "***********************************************************************");
	printWrite(caferrorLog, 1, "Global error prediction complete. Beginning error minimization for individual species.");

	moreErrSpec = [];

	for curSpec in errSpec:
		printWrite(caferrorLog, 1, "--------------------------");
		printWrite(caferrorLog, 1, "**Attempting to Minimize " + curSpec + " by adding more error");

		if Mode == 0 or Mode == 2:
			specMinScore = minScore;
		if Mode == 1 or Mode == 3:
			specMinScore = overallMinScore;
		oldspecMin = minScore;

		zcheck = 0;
		minimized = 0;
		currentError = errMin + (0.1 * errMin);
		#Some initializations.

		if Mode == 1 or Mode == 3:
			specMinErr = errMin;

		if Mode == 0 or Mode == 2:
			for eachspec in mainSpecDict:
				mainSpecDict[eachspec] = errMin;

		spec_to_min[0] = curSpec;

		while minimized == 0:
		#For each species, 10% of the global min error will be iteratively added in this loop until the score ceases to decrease OR the error has reached 1.0. At which point minimized will be set
		#to 1 and the loop will exit.

			mainSpecDict[curSpec] = currentError;

			if mainSpecDict[curSpec] <= 1.0:
				while initCheck == "bad" or initCheck == "":
					numruns = cafeRun(FamFile, Tree, LamStruct, CafePath, tmpDir, caferrorLog, mainSpecDict, spec_to_min, rmlog, FilterFlag, numruns, vOpt);
					currentScore, currentLam, initCheck = getScore(currentError, spec_to_min, tmpDir, 0, caferrorLog);
					if initCheck == "bad":
						rmlog = 1;
				initCheck = "";
				rmlog = 0;

			else:
				minimized = 1;
				continue;

			if float(currentScore) > float(specMinScore):
				check = "higher than";
				check2 = "will not be further minimized.";

				if finalSpecDict[curSpec] == 0.0:
					finalSpecDict[curSpec] = errMin;

				oldspecMin = specMinScore;
				minimized = 1;

			elif float(currentScore) == float(specMinScore):
				check = "the same as";
				check2 = "will not be further minimized.";
				oldspecMin = specMinScore;
				minimized = 1;

			elif float(currentScore) < float(specMinScore):
				check = "lower than";
				check2 = "will continue to be minimized...";

				if finalSpecDict[curSpec] == 0.0:
					finalSpecDict[curSpec] = errMin;

				oldspecMin = specMinScore;
				specMinScore = currentScore;

				if zcheck == 1:
					minimized = 1;

				finalSpecDict[curSpec] = currentError;

				specMinErr = currentError;

### 				if Mode == 1 or Mode == 3:
###					overallMinScore = currentScore;
###				Toggle this line to either keep the background constant (commented) or update it each time a species is minimized.

				currentError = currentError + (0.1 * errMin);

				if currentError >= 1:
					currentError = 1.0
					zcheck = 1;

				if curSpec not in moreErrSpec:
					moreErrSpec.append(curSpec);

		if Mode == 1 or Mode == 3:
			mainSpecDict[curSpec] = specMinErr;
			
	for curSpec in errSpec:
	#This block iteratively subtracts 10% of the global min error to each species' error model to attempt to find out which species contain LESS error than
	#the global min, and to get a general guess of what that error might be.
		if curSpec not in moreErrSpec:
			printWrite(caferrorLog, 1, "--------------------------");
			printWrite(caferrorLog, 1, "**Attempting to Minimize " + curSpec + " by adding less error");

			if Mode == 0 or Mode == 2:
				specMinScore = minScore;
			if Mode == 1 or Mode == 3:
				specMinScore = overallMinScore;
			oldspecMin = minScore;

			zcheck = 0;
			minimized = 0;
			currentError = errMin - (0.1 * errMin);

			if Mode == 1 or Mode == 3:
				specMinErr = errMin;

			if Mode == 0 or Mode == 2:
				for eachspec in mainSpecDict:
					mainSpecDict[eachspec] = errMin;

			spec_to_min[0] = curSpec;

			while minimized == 0:
			# For each species, 10% of the global min error will be iteratively subtracted in this loop until the score ceases to decrease OR the error has reached 0.0.
			# At which point minimized will be set to 1 and the loop will exit.
				mainSpecDict[curSpec] = currentError;
				
				if mainSpecDict[curSpec] >= 0.0:
					while initCheck == "bad" or initCheck == "":
						numruns = cafeRun(FamFile, Tree, LamStruct, CafePath, tmpDir, caferrorLog, mainSpecDict, spec_to_min, rmlog, FilterFlag, numruns, vOpt);
						currentScore, currentLam, initCheck = getScore(currentError, spec_to_min, tmpDir, 0, caferrorLog);
						if initCheck == "bad":
							rmlog = 1;
					initCheck = "";
					rmlog = 0;
				else:
					minimized = 1;
					continue;
				
				if float(currentScore) > float(specMinScore):
					check = "higher than";
					check2 = "will not be further minimized.";

					if finalSpecDict[curSpec] == 0.0:
						finalSpecDict[curSpec] = errMin;

					oldspecMin = specMinScore;
					minimized = 1;

				elif float(currentScore) == float(specMinScore):
					check = "the same as";
					check2 = "will not be further minimized.";

					if finalSpecDict[curSpec] == 0.0:
						finalSpecDict[curSpec] = errMin;

					oldspecMin = specMinScore;
					minimized = 1;

				elif float(currentScore) < float(specMinScore):
					check = "lower than";
					check2 = "will continue to be minimized...";
					oldspecMin = specMinScore;
					specMinScore = currentScore;

					if zcheck == 1:
						minimized = 1;

					finalSpecDict[curSpec] = currentError;

## #					if Mode == 1 or Mode == 3:
					specMinErr = currentError;
###					Toggle this line to either keep the background constant (commented) or update it each time a species is minimized.

		
					currentError = currentError - (0.1 * errMin);
	
					if currentError <= 0:
						currentError = 0.0;
						zcheck = 1;
					
			if Mode == 1 or Mode == 3:
				mainSpecDict[curSpec] = specMinErr;

# Individual species error estimation ends here if -s 1
####################################
	#####
	# The final CAFE run which generates a report with the minimized error.
	printWrite(caferrorLog, 0, "\n***********************************************************************");

	spec_to_min[0] = "final";

	while initCheck == "bad" or initCheck == "":
		numrnuns = cafeRun(FamFile, Tree, LamStruct, CafePath, tmpDir, caferrorLog, finalSpecDict, spec_to_min, rmlog, FilterFlag, numruns, vOpt);
		finalScore, finalLam, initCheck = getScore("final", spec_to_min, tmpDir, 0, caferrorLog);
		if initCheck == "bad":
			rmlog = 1;
	initCheck = "";
	rmlog = 0;
	#####
####################################
if indSpecMin == 1:
# Output block if -s 1
	printWrite(caferrorLog, 1, "# =======================================================================");
	printWrite(caferrorLog, 1, "# Final error estimates by species:");
	
	for key in finalSpecDict:
		printWrite(caferrorLog, 1, "# " + key, str(finalSpecDict[key]), pad);

	printWrite(caferrorLog, 1, "# Score with individual errors:", str(finalScore), pad);
	printWrite(caferrorLog, 1, "# Lambda with individual errors:", str(finalLam), pad);

else:
# Output block if -s 0
	spec_to_min[0] = "final";

	for eachspec in mainSpecDict:
		mainSpecDict[eachspec] = errMin;

	while initCheck == "bad" or initCheck == "":
		numrnuns = cafeRun(FamFile, Tree, LamStruct, CafePath, tmpDir, caferrorLog, mainSpecDict, spec_to_min, rmlog, FilterFlag, numruns, vOpt);
		finalScore, finalLam, initCheck = getScore("final", spec_to_min, tmpDir, 0, caferrorLog);
		if initCheck == "bad":
			rmlog = 1;
	initCheck = "";
	rmlog = 0;

####################################
# Main output block
printWrite(caferrorLog, 1, "# =======================================================================");
if firstRun == 1:
	printWrite(caferrorLog, 1, "# ************************************");
	printWrite(caferrorLog, 1, "# Score with no errormodel:", str(preScore), pad);
	printWrite(caferrorLog, 1, "# Lambda with no errormodel:", str(preLam), pad);
printWrite(caferrorLog, 1, "# ************************************");
printWrite(caferrorLog, 1, "# Global Error Estimation:", str(errMin), pad);
printWrite(caferrorLog, 1, "# Score with global errormodel:", str(minScore), pad);
printWrite(caferrorLog, 1, "# Lambda with global errormodel:", str(minLam), pad);
printWrite(caferrorLog, 1, "# ************************************");
printWrite(caferrorLog, 1, "# =======================================================================");

endsec = time.time();
runtime = (float(endsec) - float(startsec)) / float(60.0);
end = datetime.datetime.now().strftime("%m.%d.%Y | %I:%M:%S");

printWrite(caferrorLog, 1, "# Caferror finished at:", str(end), pad);
printWrite(caferrorLog, 1, "# Runtime:", str(runtime) + " minutes", pad);
#############################################################################
#End Main Block
#############################################################################





