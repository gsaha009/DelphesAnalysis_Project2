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
    
    args = parser.parse_args()

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
        'Bkg_TTJets_DiLep'     :     [109.84, '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2017/TTJets_DiLept'],
        'Bkg_TTWJetsToLNu'     :     [0.254,  '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2017/TTWJetsToLNu'],
        'Bkg_TTZJetsToLL'      :     [0.240,  '/home/gsaha/Data/TTZJetsToLL'],
        'Bkg_WZTo3LNu_012Jets' :     [2.273,  '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/WZTo3LNu_012J'],
        'Bkg_WZTo2L2Q_012Jets' :     [4.504,  '/home/gsaha/Data/WZTo2L2Q_012J_14TeV'],
        'Bkg_ZZTo4L_012Jets'   :     [0.187,  '/home/gsaha/Data/ZZTo4L_14TeV'],
        'Bkg_WWW'              :     [0.2362, '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/WWW'],
        'Bkg_WWZ'              :     [0.1889, '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/WWZ'],
        'Bkg_WZZ'              :     [0.06376,'/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/WZZ'],
        'Bkg_ZZZ'              :     [0.0158, '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/ZZZ']
    }

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
    delphesCmdList = []
    run_index = 0

    for key, valList in SampleDict.items():
        print(f'process : {key}')
        #handle = key+'_'+'MX_'+mx+'_MH_'+mh+'_'+chan+'_'+tag if 'Sig_' in key else key+'_'+chan+'_'+tag
        handle = f'{key}_MX_{mx}_MH_{mh}_{chan}_{tag}' if 'Sig_' in key else f'{key}_{chan}_{tag}'
        jobFile = os.path.join(jobcardpath, handle+'.job') 
        shutil.copy(os.path.join(os.getcwd(),'sample.job.tmpl'), jobFile)
        os.chmod(jobFile, 0o755)
        mvaInFileName  = os.path.join(joboutpath, handle+'_mvaIn.root')
        mvaOutFileName = os.path.join(joboutpath, handle+'_mvaOut.root')
        histFile       = os.path.join(joboutpath, handle+'_hist.root')
        runlogfile     = os.path.join(joboutpath, handle+'.runlog')
        
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
            for item in files:
                outf.write('inputFile '+item+'\n')
            outf.write('END')
        delphesCmdList.append(f'{EXE} {jobFile} > {runlogfile} 2>&1 &')

    with open(cmdfile, 'w') as cmdf:
        #cmdf.write('#!/bin/sh \n\n')
        for cmd in delphesCmdList:
            cmdf.write(cmd+'\n')
    

if __name__ == "__main__":
        main()
