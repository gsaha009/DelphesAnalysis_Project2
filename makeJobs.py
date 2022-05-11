import sys
import os
import glob
import argparse
import shutil
import fileinput
import ROOT
from prettytable import PrettyTable
import yaml
from stat import S_IREAD, S_IRGRP, S_IROTH
from collections import defaultdict
from subprocess import Popen, PIPE
import random

def runShellCmd(cmdList):
    proc = Popen(cmdList, stdout=PIPE)
    while True:
        output = proc.stdout.readline()
        if proc.poll() is not None:
            break
        if output:
            print(output.strip().decode("utf-8"))
    rc = proc.poll()

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)
                                            
def main():
    parser = argparse.ArgumentParser(description='JobCardMaker')
    parser.add_argument('-p',  '--proc',        type=str, action='store',      required=False, default='TTbarFCNC', help='name of the process {to make a folder tag}')
    parser.add_argument('-ap', '--anapath',     type=str, action='store',      required=False,                      help='Main analysis directory')
    parser.add_argument('-mx', '--MX',          type=str, action='store',      required=True,                       help='Mass of X')
    parser.add_argument('-mH', '--MH',          type=str, action='store',      required=True,                       help='Mass of H')
    parser.add_argument('-c',  '--channel',     type=str, action='store',      required=True,                       help='DL | SL')
    parser.add_argument('-t',  '--tag',         type=str, action='store',      required=True,                       help='CUT | BDTD | BDTG')
    parser.add_argument('-s',  '--skim',                  action='store_true', required=False,                      help='need to skim?')
    parser.add_argument('-so', '--SignalOnly',            action='store_true', required=False,                      help='jobs only for signal samples')
    parser.add_argument('-bo', '--BkgOnly',               action='store_true', required=False,                      help='jobs only for bkg samples')
    parser.add_argument('-mj', '--makejobs',              action='store_true', required=False,                      help='make job cards only')
    parser.add_argument('-pr', '--postprocess',           action='store_true', required=False,                      help='hadd output root files')
    parser.add_argument('-fr', '--force',                 action='store_true', required=False,                      help='force hadd')

    args = parser.parse_args()

    postproc = args.postprocess
    makejobs = args.makejobs

    if makejobs:
        print('Start making jobcards ...')
    elif postproc:
        print('Start doing hadd ...')
    else:
        raise RuntimeError('mention either makejobs or postprocess ... !')

    if args.SignalOnly and args.BkgOnly:
        raise RuntimeError ('Do not mention both of -so and -bo. Use none or either one of them. none will take care of the bkg samples with the signals')

    mx      = args.MX
    mh      = args.MH
    chan    = args.channel
    tag     = args.tag

    subtag  = 'Wlnu' if chan == 'DL' else 'Wjj'
    SampleDict = {
        'Sig_Xmuta_0p001_Xct_0p001_Xuc_0p005'  : [0.001, f'/store/user/gsaha/Data/FCNC_TTbar_Signals/S3_Xmuta_{subtag}_MX_{mx}_MH_{mh}_gXmuta_0p001_gXct_0p001_gXuc_0p005_DelphesRootFiles'],
        'Sig_Xmuta_0p001_Xct_0p01_Xuc_0p005'   : [0.001, f'/store/user/gsaha/Data/FCNC_TTbar_Signals/S3_Xmuta_{subtag}_MX_{mx}_MH_{mh}_gXmuta_0p001_gXct_0p01_gXuc_0p005_DelphesRootFiles'],
        'Sig_Xmuta_0p01_Xct_0p001_Xuc_0p005'   : [0.001, f'/store/user/gsaha/Data/FCNC_TTbar_Signals/S3_Xmuta_{subtag}_MX_{mx}_MH_{mh}_gXmuta_0p01_gXct_0p001_gXuc_0p005_DelphesRootFiles'],
        'Sig_Xmuta_0p01_Xct_0p01_Xuc_0p005'    : [0.001, f'/store/user/gsaha/Data/FCNC_TTbar_Signals/S3_Xmuta_{subtag}_MX_{mx}_MH_{mh}_gXmuta_0p01_gXct_0p01_gXuc_0p005_DelphesRootFiles'],
        'Sig_Xmuta_0p003_Xct_0p003_Xuc_0p005'  : [0.001, f'/store/user/gsaha/Data/FCNC_TTbar_Signals/S3_Xmuta_{subtag}_MX_{mx}_MH_{mh}_gXmuta_0p003_gXct_0p003_gXuc_0p005_DelphesRootFiles'],
        'Sig_Xmuta_0p003_Xct_0p007_Xuc_0p005'  : [0.001, f'/store/user/gsaha/Data/FCNC_TTbar_Signals/S3_Xmuta_{subtag}_MX_{mx}_MH_{mh}_gXmuta_0p003_gXct_0p007_gXuc_0p005_DelphesRootFiles'],
        'Sig_Xmuta_0p005_Xct_0p005_Xuc_0p005'  : [0.001, f'/store/user/gsaha/Data/FCNC_TTbar_Signals/S3_Xmuta_{subtag}_MX_{mx}_MH_{mh}_gXmuta_0p005_gXct_0p005_gXuc_0p005_DelphesRootFiles'],
        'Sig_Xmuta_0p007_Xct_0p003_Xuc_0p005'  : [0.001, f'/store/user/gsaha/Data/FCNC_TTbar_Signals/S3_Xmuta_{subtag}_MX_{mx}_MH_{mh}_gXmuta_0p007_gXct_0p003_gXuc_0p005_DelphesRootFiles'],
        'Sig_Xmuta_0p007_Xct_0p007_Xuc_0p005'  : [0.001, f'/store/user/gsaha/Data/FCNC_TTbar_Signals/S3_Xmuta_{subtag}_MX_{mx}_MH_{mh}_gXmuta_0p007_gXct_0p007_gXuc_0p005_DelphesRootFiles'],
        ##'Bkg_TTJets_DiLep_Set_Old'             : [109.84, '/store/user/gsaha/Data/Backgrounds/TTJets_DiLept'],
        'Bkg_TTJets_DiLep_Set_New'             : [107.65, '/store/user/gsaha/Data/Backgrounds/TTbarDL_elmuta_012Jets_LO_MLM_DelphesRootFiles'],
        'Bkg_TTJets_SingleLep'                 : [437.14, '/store/user/gsaha/Data/Backgrounds/TTbarSL_elmuta_012Jets_LO_MLM_DelphesRootFiles'],
        'Bkg_TTWJetsToLNu'                     : [0.254,  '/store/user/gsaha/Data/Backgrounds/TTWJetsToLNu'],
        'Bkg_TTWToQQ'                          : [0.103,  '/store/user/gsaha/Data/Backgrounds/TTWnojetsToQQ_LO_DelphesRootFiles'],
        'Bkg_TTZJetsToLL'                      : [0.240,  '/store/user/gsaha/Data/Backgrounds/TTZJetsToLL'],
        'Bkg_WZTo3LNu_012Jets'                 : [2.273,  '/store/user/gsaha/Data/Backgrounds/WZTo3LNu_012J'],
        'Bkg_WZTo2L2Q_012Jets'                 : [4.504,  '/store/user/gsaha/Data/Backgrounds/WZTo2L2Q_012J_14TeV'],
        'Bkg_ZZTo4L_012Jets'                   : [0.187,  '/store/user/gsaha/Data/Backgrounds/ZZTo4L_14TeV'],
        'Bkg_WWW'                              : [0.2362, '/store/user/gsaha/Data/Backgrounds/WWW'],
        'Bkg_WWZ'                              : [0.1889, '/store/user/gsaha/Data/Backgrounds/WWZ'],
        'Bkg_WZZ'                              : [0.06376,'/store/user/gsaha/Data/Backgrounds/WZZ'],
        'Bkg_ZZZ'                              : [0.0158, '/store/user/gsaha/Data/Backgrounds/ZZZ'],
        'Bkg_bbtautau'                         : [0.114,  '/store/user/gsaha/Data/Backgrounds/bbtautau_QCD_QED_LO_DelphesRootFiles'],
        'Bkg_TTHnojetsToTauTau'                : [0.006,  '/store/user/gsaha/Data/Backgrounds/TTHnojetsToTauTau_LO_DelphesRootFiles'],
        'Bkg_TTZnojetsToQQ'                    : [0.206,  '/store/user/gsaha/Data/Backgrounds/TTZnojetsToQQ_LO_DelphesRootFiles']
    }
    nodes = ['compute-0-1', 'compute-0-2', 'compute-0-3', 'compute-0-4', 'compute-0-6', 'compute-0-7', 'compute-0-8', 'compute-0-9', 'compute-0-10']
    #'compute-0-5', 'compute-0-15', 'compute-0-16', 'compute-0-17', 'compute-0-18', 'compute-0-19','compute-0-20','compute-0-21']

    analysispath = os.getcwd() if args.anapath == None else args.anapath
    jobcardpath  = os.path.join(analysispath, f'JobCards_MX_{mx}') 
    joboutpath   = os.path.join(analysispath, f'JobOutput_MX_{mx}')
    EXE          = os.path.join(analysispath, 'FCNC.exe') 
    
    if not os.path.isdir(jobcardpath):
        os.mkdir(jobcardpath)
    else:
        print(f'{jobcardpath} exists !!!')

    if not os.path.isdir(joboutpath):
        os.mkdir(joboutpath)
    else:
        print(f'{joboutpath} exists !!!')

    #cmdfile = os.path.join(jobcardpath, 'DelphesCmds'+'_MX_'+mx+'_MH_'+mh+'_'+chan+'_'+tag+'.log')
    cmdfile = os.path.join(jobcardpath, f'DelphesCmds_MX_{mx}_MH_{mh}_{chan}_{tag}.log')
    sgefile = os.path.join(jobcardpath, f'DelphesCmds_MX_{mx}_MH_{mh}_{chan}_{tag}.sh')
    sgerefile = os.path.join(jobcardpath, f'DelphesCmds_MX_{mx}_MH_{mh}_{chan}_{tag}_Resubmit.sh')
    delphesCmdList = []
    sgeCmdList = []
    sgeResubCmdList = []
    run_index = 0
    histfilestoadddict = defaultdict(list)
    treefilestoadddict = defaultdict(list)
    haddfile = os.path.join(jobcardpath, f'DelphesCmds_MX_{mx}_MH_{mh}_{chan}_{tag}_hadd.sh')
    haddf = open(haddfile, 'w')
    haddf.write('#!/bin/sh \n\n')

    for key, valList in SampleDict.items():
        print(f'process : {key}')
        handle = f'{key}_MX_{mx}_MH_{mh}_{chan}_{tag}' if 'Sig_' in key else f'{key}_{chan}_{tag}'
        
        xsec           = valList[0]
        filepath       = valList[1]
        files          = glob.glob(os.path.join(filepath, '*.root'))

        nEntries = 0
        for file in files:
            print(f'\t infile : {file}')
            tfile = ROOT.TFile(file, 'read')
            if not tfile or tfile.IsZombie():
                print(f'WARNING : {file} is zombie !!!')
                continue
            treeList = [key.GetName() for key in tfile.GetListOfKeys()]
            if len(treeList) < 2 or treeList[1] != 'Delphes':
                print('TTree :Delphes: not found !!')
                continue
            ttree = tfile.Get('Delphes')
            nEntries += ttree.GetEntries()

        print(f'\t ... Cross section : {xsec} pb, nEvents : {nEntries}')    

        BDT_WEIGHT1 = f'/home/gsaha/Work/DelphesML/Project2/BDT/IC_Model_ChiMass_{mx}_v2/weights/TMVAClassification_{args.tag}.weights.xml'

        filesPerJob = 5
        infileListPerJob = [files[i:i+filesPerJob] for i in range(0, len(files), filesPerJob)]

        #tempdict = defaultdict(list)
        
        for i, filelist in enumerate(infileListPerJob):
            mvaInFileName  = os.path.join(joboutpath, f'{handle}_{i}_mvaIn.root')
            mvaOutFileName = os.path.join(joboutpath, f'{handle}_{i}_mvaOut.root')
            histFile       = os.path.join(joboutpath, f'{handle}_{i}_hist.root')
            runlogfile     = os.path.join(joboutpath, f'{handle}_{i}.runlog')
            _outfile       = os.path.join(joboutpath, f'{handle}_{i}.out')
            _errfile       = os.path.join(joboutpath, f'{handle}_{i}.err')
            jobFile = os.path.join(jobcardpath, f'{handle}_{i}.job')
            sgeFile = os.path.join(jobcardpath, f'{handle}_{i}.sh')
                
            if makejobs:
                shutil.copy(os.path.join(os.getcwd(),'sample.job.tmpl'), jobFile)
                os.chmod(jobFile, 0o755)
                shutil.copy(os.path.join(os.getcwd(),'sample.sh.tmpl'), sgeFile)
                os.chmod(sgeFile, 0o755)

                # Making job file
                replaceAll(jobFile, 'MVAnetwork BDTtype', 'MVAnetwork '+args.tag)
                replaceAll(jobFile, 'MVAxmlFile BDT_WEIGHT1', 'MVAxmlFile '+BDT_WEIGHT1)
                replaceAll(jobFile, 'mvaInputFile mvaIn_SAMPLE', 'mvaInputFile '+mvaInFileName)
                replaceAll(jobFile, 'mvaOutputFile mvaOut_SAMPLE', 'mvaOutputFile '+mvaOutFileName)
                replaceAll(jobFile, 'lumiWtList xsec=xsec_SAMPLE intLumi=lumi_SAMPLE nevents=500000', 'lumiWtList xsec='+str(xsec)+' '+'intLumi='+str(300000)+' nevents='+str(nEntries))
                replaceAll(jobFile, 'histFile  hist_SAMPLE', 'histFile  '+histFile)
            
                if args.tag == "CUT":
                    replaceAll(jobFile, 'cutbased val', 'cutbased 1')
                    replaceAll(jobFile, 'mvabased val', 'mvabased 0')
                    replaceAll(jobFile, 'createMVATree val', 'createMVATree 1') if args.skim else replaceAll(jobFile, 'createMVATree val', 'createMVATree 0')
                    replaceAll(jobFile, 'readMVA val', 'readMVA 0')
                elif args.tag == "BDTD" or args.tag =="BDTG":
                    replaceAll(jobFile, 'cutbased val', 'cutbased 0')
                    replaceAll(jobFile, 'mvabased val', 'mvabased 1')
                    replaceAll(jobFile, 'createMVATree val', 'createMVATree 1') if args.skim else replaceAll(jobFile, 'createMVATree val', 'createMVATree 0')
                    replaceAll(jobFile, 'readMVA val', 'readMVA 1') if not args.skim else replaceAll(jobFile, 'readMVA val', 'readMVA 0') 
                else:
                    replaceAll(jobFile, 'createMVATree val', 'createMVATree 0')
                    replaceAll(jobFile, 'readMVA val', 'readMVA 0')            
                    raise RuntimeError("Mention BDT or CUT")

                replaceAll(jobFile, 'lowXmass val', 'lowXmass 1')
                replaceAll(jobFile, 'highXmass val', 'highXmass 0')
                replaceAll(jobFile, 'bdtThreshold score', 'bdtThreshold -0.03')
        
                if chan == 'DL':
                    replaceAll(jobFile, 'isDL val', 'isDL 1')
                    replaceAll(jobFile, 'isSL val', 'isSL 0')
                elif chan == 'SL':
                    replaceAll(jobFile, 'isDL val', 'isDL 0')
                    replaceAll(jobFile, 'isSL val', 'isSL 1')

                with open(jobFile, 'a') as outf:
                    for item in filelist:
                        outf.write('inputFile '+item+'\n')
                    outf.write('END')

                delphesCmdList.append(f'{EXE} {jobFile} > {runlogfile} 2>&1 &')

                # Making sge script for each job
                replaceAll(sgeFile, '__name', f'{handle}_{i}')
                replaceAll(sgeFile, '__node', f'{random.choice(nodes)}')
                replaceAll(sgeFile, '__jobcardpath', f'{jobcardpath}')
                replaceAll(sgeFile, '__jobout.output', f'{_outfile}')
                replaceAll(sgeFile, '__joberr.output', f'{_errfile}')
                replaceAll(sgeFile, '__jobcardbasename', f'{os.path.basename(jobFile)}')

                sgeCmdList.append(['qsub', sgeFile])


            elif postproc:
                issue = False
                if not os.path.exists(histFile):
                    print(f'{histFile} found missing. Plese reproduce the file')
                    issue = True
                else :
                    histf = ROOT.TFile(histFile, 'r')
                    if histf.IsZombie():
                        print(f'{histFile} is zombie. Plese reproduce the file')
                        issue = True
                    else:
                        histfilestoadddict[key].append(histFile)

                if args.skim:
                    if not os.path.exists(mvaInFileName):
                        print(f'{mvaInFileName} found missing. Plese reproduce the file')
                        issue = True
                    else:
                        treef = ROOT.TFile(mvaInFileName, 'r')
                        if histf.IsZombie():
                            print(f'{mvaInFileName} is zombie. Plese reproduce the file')
                            issue = True
                        else:
                            treefilestoadddict[key].append(mvaInFileName)
                
                if issue : sgeResubCmdList.append(['qsub', sgeFile])

            else:
                raise RuntimeError('mention makejobs or postproc')



    if makejobs:
        with open(cmdfile, 'w') as cmdf:
            #cmdf.write('#!/bin/sh \n\n')
            for cmd in delphesCmdList:
                cmdf.write(cmd+'\n')
        with open(sgefile, 'w') as sgef:
            sgef.write('#!/bin/sh \n\n') 
            for cmd in sgeCmdList:
                sgef.write(' '.join(cmd)+'\n')
    
    
    elif postproc:
        if len(sgeResubCmdList) == 0 or args.force:
            print('No failed jobs. Starting to postprocess')
            hadddir = os.path.join(joboutpath, 'results')
            if not os.path.isdir(hadddir):
                os.mkdir(hadddir)
            else:
                print(f'{hadddir} dir exists ...')

            histhaddcmds = []
            treehaddcmds = []
            for key, vallist in histfilestoadddict.items():
                handle = f'{key}_MX_{mx}_MH_{mh}_{chan}_{tag}' if 'Sig_' in key else f'{key}_{chan}_{tag}'
                histfiletoprod = f'{handle}_hist.root'
                histhaddcmds = ['hadd',os.path.join(hadddir,histfiletoprod)] + vallist
                #print('.............. '+ histhaddcmds)
                haddf.write(f'# adding hist files : {handle} ---> nFiles : {len(vallist)}\n')
                haddf.write(' '.join(histhaddcmds)+'\n\n')
                #runShellCmd(histhaddcmds)
                
            if args.skim:
                print('Please mention -s if you want to add the tree files')
                for key, vallist in treefilestoadddict.items():
                    handle = f'{key}_MX_{mx}_MH_{mh}_{chan}_{tag}' if 'Sig_' in key else f'{key}_{chan}_{tag}'
                    treefiletoprod = f'{handle}_mvaIn.root'
                    treehaddcmds = ['hadd',os.path.join(hadddir,treefiletoprod)] + vallist
                    #print('.............. '+ treehaddcmds)
                    #haddf.write('#!/bin/sh \n\n')
                    haddf.write(f'# adding tree files : {handle} ---> nFiles : {len(vallist)}\n')
                    haddf.write(' '.join(treehaddcmds)+'\n\n')
                    #runShellCmd(treehaddcmds)

        else: # there are failed jobs i.e. files are missing
            with open(sgerefile, 'w') as sgefr:
                sgefr.write('#!/bin/sh \n\n') 
                for cmd in sgeResubCmdList:
                    sgefr.write(' '.join(cmd)+'\n')
            print('Resubmit these jobs and do the postprocess again')
    
    else:
        raise RuntimeError('Mention makejobs or postprocess ... !')


if __name__ == "__main__":
        main()
