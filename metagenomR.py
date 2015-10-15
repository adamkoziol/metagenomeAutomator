__author__ = 'akoziol'

import os, errno, re, shutil, subprocess, json, sys, time, gzip
from glob import glob
from argparse import ArgumentParser
from multiprocessing import Pool
from collections import defaultdict

#Parser for arguments
parser = ArgumentParser(description='Prep Illumina fastq metagenome files to be processed by OneCodex')
parser.add_argument('-v', '--version', action='version', version='%(prog)s')
parser.add_argument('-p', '--path', required=True, help='Specify path')
parser.add_argument('-u', '--upload', required=False, default=False, help='Upload metagenome files following preparation steps?')

# Get the arguments into a list
args = vars(parser.parse_args())

# Define variables from the arguments - there may be a more streamlined way to do this
path = args['path']
upload = args['upload']
apikey = '17d05d3bfde945fea1f3d1a1db5a0b1f'

# Start time
start = time.time()

# Welcome message
print("Welcome to the CFIA Metagenome preparation pipeline.")


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

# Initialise the dictionary responsible for storing the report data
metadataFile = defaultdict(make_dict)


def commandR(sampleNames, path):
    """Opens the *name*_metadataCollection.json file and extracts any commands performed from previous iterations of the
    pipeline on data in the supplied path"""
    # Initialise the command dictionary
    performedCommands = defaultdict(make_dict)
    # Open the *name*_metadataCollection.json file for each sample
    for name in sampleNames:
        countSize = 0
        if os.path.isfile("%s/%s/%s_metadataCollection.json" % (path, name, name)):
            countSize = os.stat("%s/%s/%s_metadataCollection.json" % (path, name, name)).st_size
            if countSize != 0:
                with open("%s/%s/%s_metadataCollection.json" % (path, name, name)) as jsonReport:
                    # Load the data
                    jsonData = json.load(jsonReport)
                    if jsonData:
                        # Find the precise command used
                        for command in jsonData["commands"]:
                            # If the command exists, and is in the right format
                            if not "N/A" in jsonData["commands"][command] and not "defaultdict" in jsonData["commands"][command]:
                                # Populate the performedCommands dictionary as appropriate
                                performedCommands[name][str(command)] = str(jsonData["commands"][command])
    # Return the dictionary
    return performedCommands


def jsonR(sampleNames, path, metadata, fileType):
    """Creates a JSON report from a supplied metadata file"""
    # Make the reports folder as required
    reportPath = "%s/reports" % path
    make_path(reportPath)
    for name in sampleNames:
            # Create the .json file for each sample
            newPath = path + "/" + name
            reportName = "%s_metadata%s.json" % (name, fileType)
            JSONreport = open("%s/%s" % (newPath, reportName), "wb")
            # Print the JSON data to file
            output = json.dumps(metadata[name], sort_keys=True, indent=4, separators=(',', ': '))
            JSONreport.write(output)
            JSONreport.close()
            # Move all the reports to a common directory
            shutil.copy("%s/%s" % (newPath, reportName), "%s/%s" % (reportPath, reportName))


def foldererPrepProcesses(sampleName, path):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    foldererPrepArgs = []
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == '__main__':
        createfoldererPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            foldererPrepArgs.append((name, path))
        # This map function allows for multi-processing
        createfoldererPool.map(folderer, foldererPrepArgs)


def folderer((name, path)):
    """Uses gzip to decompress .gz.fastq files"""
    # Get any .gz file starting with the strain name into a list
    gzFiles = glob("%s/%s/*.gz" % (path, name))
    # If there are any .gz files...
    if gzFiles:
        for gzFile in gzFiles:
            # Gzip each file
            gzipCommand = "gzip -d --force %s" % gzFile
            subprocess.call(gzipCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    sys.stdout.write('.')


def fileR(list):
    """Helper script that creates a set of the stain names created by stripping off parts of the filename.
    Hopefully handles different naming conventions (e.g. 2015-SEQ-001_S1_L001_R1_001.fastq(.gz),
    2015-SEQ-001_R1_001.fastq.gz, 2015-SEQ-001_R1.fastq.gz, 2015-SEQ-001_1.fastq.gz, and 2015-SEQ-001_1.fastq.gz
    all become 2015-SEQ-001)"""
    # Initialise the set
    fileSet = set()
    for seqFile in list:
        # Search for the conventional motifs present following strain names
        # _S\d+_L001_R\d_001.fastq(.gz) is a typical unprocessed Illumina fastq file
        if re.search("_S\d+_L001", seqFile):
            fileSet.add(re.split("_S\d+_L001", seqFile)[0])
        # Files with _R\d_001.fastq(.gz) are created in the SPAdes assembly pipeline
        elif re.search("_R\d_001", seqFile):
            fileSet.add(re.split("_R\d_001", seqFile)[0])
        # _R\d.fastq(.gz) represents a simple naming scheme for paired end reads
        elif re.search("R\d.fastq", seqFile):
            fileSet.add(re.split("_R\d.fastq", seqFile)[0])
        # _\d.fastq(.gz) is always possible
        elif re.search("[-_]\d.fastq", seqFile):
            fileSet.add(re.split("[-_]\d.fastq", seqFile)[0])
        # .fastq is the last option
        else:
            fileSet.add(re.split(".fastq", seqFile)[0])
        sys.stdout.write('.')
    return fileSet


def nameExtractor():
    """Creates a set of the names of the files to be used in the analyses"""
    fileSet = set()
    # Create lists of the .gz, the .fastq, and the folders in the path
    gzChecker = glob("*.gz")
    fastqChecker = glob("*.fastq")
    folderChecker = glob("*/")
    # For each list, ensure that the list exists...
    if gzChecker:
        # Extract the unique names using the helper function fileR...
        fileList = fileR(gzChecker)
        # Add the results to fileset
        for seqFile in fileList:
            fileSet.add(seqFile)
    if fastqChecker:
        fileList = fileR(fastqChecker)
        for seqFile in fileList:
            fileSet.add(seqFile)
    if folderChecker:
        for seqFolder in folderChecker:
            # Exclude the 'reports' folder from the analysis
            if not "reports" in seqFolder and not "results" in seqFolder:
                seqName = seqFolder.rsplit("/")[0]
                fileSet.add(seqName)
    # Create a more easily parseable list from the set
    fileList = list(fileSet)
    # Return it
    return fileList


def seqMovR(path, seqNames):
    """Creates appropriately named folders, and moves sequence files to appropriate folders"""
    for seqName in seqNames:
        # Make folders as required
        make_path("%s/%s" % (path, seqName))
        # Search the path for any file or folder that contains the seqName
        filecheck = [f for f in os.listdir(path) if re.search("%s" % seqName, f)]
        for seqFile in filecheck:
            # Move files, ignore folders
            if os.path.isfile(seqFile):
                sys.stdout.write('.')
                shutil.move(seqFile, "%s/%s/%s" % (path, seqName, seqFile))
        sys.stdout.write('.')


def trimmomaticPrepProcesses(sampleName, path, metadata, commands):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    trimmomaticPrepArgs = []
    output = {}
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == '__main__':
        createtrimmomaticPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            trimmomaticPrepArgs.append((name, path, metadata, commands))
        # This map function allows for multi-processing
        output = createtrimmomaticPool.map(trimmomatic, trimmomaticPrepArgs)
    return output


def trimmomatic((name, path, metadata, commands)):
    """Runs trimmomatic to trim .fastq files on phred sequence quality. Populates metadata files as required"""
    # Get the .fastq files in the sample folder
    seqFiles = sorted(glob("%s/%s/*.fastq" % (path, name)))
    # create a handy variable for use in the trimmomatic system call
    newPath = "%s/%s/%s" % (path, name, name)
    # Ensure that trimmomatic has not already been run on these files
    if not commands[name]["trimmomaticCall"]:
        # Treat paired-end samples slightly differently than single-end samples
        if len(seqFiles) == 2:
            # Prepare the trimmomatic call
            trimmomaticCall = "trimmomatic-0.30.jar PE -threads 24 -phred33 -trimlog %s.log " \
                "%s %s %s.paired1.fq %s.unpaired1.fq %s.paired2.fq %s.unpaired2.fq " \
                "ILLUMINACLIP:/home/blais/Bioinformatics/Trimmomatic-0.30/adapters/TruSeq3-PE.fa:2:30:10 " \
                "LEADING:10 TRAILING:10 SLIDINGWINDOW:4:30 MINLEN:36" \
                % (newPath, seqFiles[0], seqFiles[1], newPath, newPath, newPath, newPath)
        else:
            trimmomaticCall = "trimmomatic-0.30.jar SE -threads 24 -phred33 -trimlog %s.log " \
                      "%s %s.trimmed.fq " \
                      "ILLUMINACLIP:/home/blais/Bioinformatics/Trimmomatic-0.30/adapters/TruSeq3-PE.fa:2:30:10 " \
                      "LEADING:10 TRAILING:10 SLIDINGWINDOW:4:30 MINLEN:36" \
                      % (newPath, seqFiles[0], newPath)
        # Run trimmomatic
        subprocess.call(trimmomaticCall, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        # Populate the metadata dictionary with the system call
        metadata[name]["commands"]["trimmomaticCall"] = trimmomaticCall
        # Remove the original .fastq files
        for seqFile in seqFiles:
            os.remove(seqFile)
    else:
        metadata[name]["commands"]["trimmomaticCall"] = commands[name]["trimmomaticCall"]
    sys.stdout.write('.')
    return metadata


def pairedEndJoinerPrepProcesses(sampleName, path, metadata, commands):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    pairedPrepArgs = []
    output = {}
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == '__main__':
        createpairedPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            pairedPrepArgs.append((name, path, metadata, commands))
        # This map function allows for multi-processing
        output = createpairedPool.map(endPairer, pairedPrepArgs)
    return output


def endPairer((name, path, metadata, commands)):
    """Joins quality-trimmed paired-end fastq files"""
    seqFiles = sorted(glob("%s/%s/*.paired*" % (path, name)))
    # create a handy variable for use in the endpairer system call
    newPath = "%s/%s/%s" % (path, name, name)
    # Ensure that join_paired_ends/py has not already been run on these files
    if not commands[name]["endPairerCall"]:
        # Treat paired-end samples slightly differently than single-end samples
        if len(seqFiles) == 2:
            # endPairerCall = "join_paired_ends.py %s %s > %s.trimmed.fq" % (seqFiles[0], seqFiles[1], newPath)
            endPairerCall = "Pear -f %s -r %s -j 24 -q 10 -o %s" % (seqFiles[0], seqFiles[1], newPath)
            subprocess.call(endPairerCall, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            # print endPairerCall
            metadata[name]["commands"]["endPairerCall"] = endPairerCall
    else:
        metadata[name]["commands"]["endPairerCall"] = commands[name]["endPairerCall"]
    # Remove the .paired\d.fq files
    for pairedFile in seqFiles:
        os.remove(pairedFile)
    sys.stdout.write('.')
    return metadata


def mergePrepProcesses(sampleName, path, metadata, commands):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    mergePrepArgs = []
    output = {}
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == '__main__':
        createmergePool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            mergePrepArgs.append((name, path, metadata, commands))
        # This map function allows for multi-processing
        output = createmergePool.map(mergR, mergePrepArgs)
    return output


def mergR((name, path, metadata, commands)):
    """Merges the joined reads and the unpaired reads"""
    make_path("%s/results" % path)
    # Use a piped cat and gzip command
    if not commands[name]["readMergeCall"]:
        readMergeCall = "cat %s/%s/*.fq %s/%s/*.fastq| gzip > %s/results/%s.fastq.gz" % (path, name, path, name, path, name)
        subprocess.call(readMergeCall, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        metadata[name]["commands"]["readMergeCall"] = readMergeCall
    else:
        metadata[name]["commands"]["readMergeCall"] = commands[name]["readMergeCall"]
    sys.stdout.write('.')
    #Remove .log, .fq and .fastq files - the folders will remain, as they will be used to get the names of the strains
    # if the pipeline needs to be run again.
    seqFiles = glob("%s/%s/*" % (path, name))
    if seqFiles:
        for seqFile in seqFiles:
            os.remove(seqFile)
    return metadata


def uploadR(sampleNames, path, metadata, commands):
    """Uploads samples to OneCodex for analysis"""
    """Not fully tested, as it seems to be impossible to delete uploaded sequences right now, and I don't want to
    fill up Cathy's OneCodex account"""
    os.chdir("%s/reports" % path)
    for seqName in sampleNames:
        if not commands["uploadCommand"]:
            curlCommand = "curl https://beta.onecodex.com/api/v0/upload -X POST -u %s: " \
                  "--form filename=@%s.fastq.gz" % (apikey, seqName)
            print curlCommand
            proc = subprocess.Popen(curlCommand, shell=True, stdout=subprocess.PIPE, stderr=open(os.devnull, 'wb'))
            stdout_value = proc.communicate()[0]
            print stdout_value
            metadata[seqName]["commands"]["uploadCommand"] = curlCommand
            metadata[seqName]["commands"]["sample_id"] = stdout_value
        else:
            metadata[seqName]["commands"]["uploadCommand"] = commands["uploadCommand"]
            metadata[seqName]["commands"]["sample_id"] = commands["sample_id"]
    os.chdir(path)
    return metadata


def filler(metadata, metadataList):
    """Properly populates the metadata dictionary - when I tried to populate the metadata dictionary within the
     multi-processed functions, it was returned as a list (likely due to how the files were returned. This essentially
     iterates through the list, and populates a dictionary appropriately"""
    # Make a copy of the metadata dictionary
    metadataCopy = metadata
    # Iterate through metadataCopy
    for item in metadataList:
        # each item in metadataCopy is a dictionary entry
        # iterate through all the dictionaries
        for name in item:
            # The way the dictionaries were created, they should have the format:
            # metadataOutput[name]["commands"]["trimmomaticCall"] = "text text text...."
            for generalCategory in item[name]:
                for specificCategory in item[name][generalCategory]:
                    # Populate the dictionary
                    if specificCategory not in metadataCopy[name][generalCategory]:
                        metadataCopy[name][generalCategory][specificCategory] = str(item[name][generalCategory][specificCategory])
    # Return the beautifully-populated dictionary
    return metadataCopy


def runMetagenomR():
    os.chdir(path)
    print "Finding sample names"
    seqNames = nameExtractor()
    print "\nMoving files to appropriate folders"
    seqMovR(path, seqNames)
    print "\nExtracting files from archive as necessary"
    foldererPrepProcesses(seqNames, path)
    commands = commandR(seqNames, path)
    print "\nPerforming trimmomatic quality trimming on samples"
    trimmomaticMetadataList = trimmomaticPrepProcesses(seqNames, path, metadataFile, commands)
    trimmomaticMetadata = filler(metadataFile, trimmomaticMetadataList)
    jsonR(seqNames, path, trimmomaticMetadata, "Collection")
    print "\nMerging paired ends and appending singletons"
    pairedMetadataList = pairedEndJoinerPrepProcesses(seqNames, path, trimmomaticMetadata, commands)
    pairedMetadata = filler(trimmomaticMetadata, pairedMetadataList)
    jsonR(seqNames, path, pairedMetadata, "Collection")
    mergeMetadataList = mergePrepProcesses(seqNames, path, pairedMetadata, commands)
    mergeMetadata = filler(pairedMetadata, mergeMetadataList)
    jsonR(seqNames, path, mergeMetadata, "Collection")
    # Upload the files to OneCodex if desired
    if upload:
        print "\nUploading sequences to OneCodex"
        uploadMetadata = uploadR(seqNames, path, mergeMetadata, commands)
        jsonR(seqNames, path, uploadMetadata, "Collection")
    # print json.dumps(commands, sort_keys=True, indent=4, separators=(',', ': '))

runMetagenomR()
print "\nElapsed Time: %.2f seconds" % (time.time() - start)

