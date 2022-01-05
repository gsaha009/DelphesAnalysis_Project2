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

def sampleInfo():
    return yaml.safe_load(
        '''
        file:
        xsec:
        nEvents:
        '''
        )

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)
                                            
def main():
    parser = argparse.ArgumentParser(description='JobCardMaker')
    parser.add_argument('--proc',        type=str, action='store',      required=False, default='TTbarFCNC', help='name of the process {XZ or HZ}')
    parser.add_argument('--channel',     type=str, action='store',      required=True,                       help='DL | SL')
    parser.add_argument('--MX',          type=str, action='store',      required=True,                       help='Mass of X')
    parser.add_argument('--tag',         type=str, action='store',      required=True,                       help='BDT | CUT | DNN ...')
    parser.add_argument('--skim',                  action='store_true', required=False,                      help='need to skim?')
    parser.add_argument('--SignalOnly',            action='store_true', required=False,                      help='jobs only for signal samples')
    parser.add_argument('--BkgOnly',               action='store_true', required=False,                      help='jobs only for bkg samples')

    args = parser.parse_args()
    
    couplList = [(0.001,0.001),(0.001,0.003),(0.001,0.005),(0.001,0.007),(0.001,0.009),(0.001,0.01),
                 (0.003,0.001),(0.003,0.003),(0.003,0.005),(0.003,0.007),(0.003,0.009),(0.003,0.01),
                 (0.005,0.001),(0.005,0.003),(0.005,0.005),(0.005,0.007),(0.005,0.009),(0.005,0.01),
                 (0.007,0.001),(0.007,0.003),(0.007,0.005),(0.007,0.007),(0.007,0.009),(0.007,0.01),
                 (0.009,0.001),(0.009,0.003),(0.009,0.005),(0.009,0.007),(0.009,0.009),(0.009,0.01),
                 (0.01,0.001),(0.01,0.003),(0.01,0.005),(0.01,0.007),(0.01,0.009),(0.01,0.01)]
    
    BkgDict = {
        'Bkg_DY1JetsToLL' :      [1188.6, '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/DY1JetsToLL'],
        #'Bkg_DY2JetsToLL' :      [523.0,  '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/DY2JetsToLL'],
        'Bkg_DY3JetsToLL' :      [169.4,  '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/DY3JetsToLL'],
        'Bkg_TTJets_DiLep' :     [109.84, '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2017/TTJets_DiLept'],
        'Bkg_TTWJetsToLNu' :     [0.254,  '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2017/TTWJetsToLNu'],
        #'Bkg_TTZJetsToLL' :      [0.662,  '/home/gsaha/sshfs_mount/sinpcms_Disk1/data/GOURAB/TTZJetsToLL'],
        'Bkg_TTZJetsToLL' :      [0.240,  '/home/gsaha/Data/TTZJetsToLL'],
        #'Bkg_WZTo3LNu_0Jets' :   [1.0,    '/home/gsaha/sshfs_mount/sinpcms_Disk1/data/GOURAB/Last7Days/BkgProductuon/WZTo3LNu_Exclusive/WZTo3LNu_0J'],
        #'Bkg_WZTo3LNu_1Jets' :   [1.0,    '/home/gsaha/sshfs_mount/sinpcms_Disk1/data/GOURAB/Last7Days/BkgProductuon/WZTo3LNu_Exclusive/WZTo3LNu_1J'],
        #'Bkg_WZTo3LNu_2Jets' :   [1.0,    '/home/gsaha/sshfs_mount/sinpcms_Disk1/data/GOURAB/Last7Days/BkgProductuon/WZTo3LNu_Exclusive/WZTo3LNu_2J'],
        #'Bkg_WZTo3LNu_3Jets' :   [1.0,    '/home/gsaha/sshfs_mount/sinpcms_Disk1/data/GOURAB/Last7Days/BkgProductuon/WZTo3LNu_Exclusive/WZTo3LNu_3J'],
        'Bkg_WZTo3LNu_012Jets' : [2.273,   '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/WZTo3LNu_012J'],
        #'Bkg_WZTo2L2Q_012Jets' : [5.6,    '/home/gsaha/sshfs_mount/sinpcms_Disk1/data/GOURAB/Last7Days/BkgProductuon/WZTo2L2Q_012J_14TeV'],
        'Bkg_WZTo2L2Q_012Jets' : [4.504,    '/home/gsaha/Data/WZTo2L2Q_012J_14TeV'],
        #'Bkg_ZZTo4L_012Jets' :   [0.23,   '/home/gsaha/sshfs_mount/sinpcms_Disk1/data/GOURAB/Last7Days/BkgProductuon/ZZTo4L_14TeV'],
        'Bkg_ZZTo4L_012Jets' :   [0.187,   '/home/gsaha/Data/ZZTo4L_14TeV'],
        'Bkg_GGF_ZZTo4L' :       [0.01482,'/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/ggfHToZZTo4L_14TeV'],
        'Bkg_VBF_ZZTo4L' :       [0.00221,'/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/vbfHToZZTo4L_14TeV'],
        'Bkg_WWW' :              [0.2362, '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/WWW'],
        'Bkg_WWZ' :              [0.1889, '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/WWZ'],
        'Bkg_WZZ' :              [0.06376,'/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/WZZ'],
        'Bkg_ZZZ' :              [0.0158, '/home/gsaha/sshfs_mount/lacie8tby/Exotic/Delphes/2016/ZZZ']
    }

    MX      = args.MX
    chan    = args.channel
    masstag = args.proc+'_ChiMass_'+str(MX)+'_'+chan

    infilepath   = '/home/gsaha/Data/SIGNALs'
    analysispath = os.getcwd()
    jobcardpath  = os.path.join(analysispath, 'JobCards_'+masstag) 
    joboutpath   = os.path.join(analysispath, 'JobOutput_'+masstag)
    EXE          = os.path.join(analysispath, 'FCNC.exe') 
    
    if not os.path.isdir(jobcardpath):
        os.mkdir(jobcardpath)
    else:
        print(f'{jobcardpath} exists !!!')

    if not os.path.isdir(joboutpath):
        os.mkdir(joboutpath)
    else:
        print(f'{joboutpath} exists !!!')

    bkgKeyList = list(BkgDict.keys())

    if args.SignalOnly:
        keyList = couplList
    elif args.BkgOnly:
        keyList = bkgKeyList
    else:
        keyList = couplList + bkgKeyList

    print(keyList)
    ptab = PrettyTable(['Process','x-sec (pb)','nEvents produced'])
    ptab.title = f'Process : {masstag}'
    logfile = os.path.join(jobcardpath, masstag+'_evInfo.log')
    yamlfile = os.path.join(jobcardpath, masstag+'_DNN.yaml')
    cmdfile = os.path.join(jobcardpath, masstag+'_delphesCmds.log')
    delphesCmdList = []
    procdict = {
        'SignalSamples': {},
        'Train_Background': {},
        'Other_Background': {}
    }
    sampledictSig = dict()
    sampledictTrainBkg = dict()
    sampledictOtherBkg = dict()
    infodictSig = dict()
    infodictTrainBkg = dict()
    infodictOtherBkg = dict()


    for key in keyList:
        isBkg   = True if 'Bkg' in key else False
        coupl   = 'Yl_'+str(key[0])+'_Yq_'+str(key[1]) if not isBkg else key
        handle  = masstag+'_'+coupl
        ptab_handle = handle if not isBkg else handle.replace(masstag+'_Bkg_','')
        infiledir = os.path.join(infilepath, os.path.join(masstag, handle+'_RootFiles')) if not isBkg else BkgDict.get(key)[1]
        bannerdir = os.path.join('/home/gsaha/Work/MadGenerate/SIGNAL_XSEC_CHECK', os.path.join(handle,'Events/run_01')) if not isBkg else BkgDict.get(key)[1]
        #bannerdir = os.path.join('/home/gsaha/Work/MadGenerate/SIGNAL_XSEC_CHECK', os.path.join(handle,'Events')) if not isBkg else BkgDict.get(key)[1]
        print(f' Sample : {handle}')
        if not os.path.isdir(infiledir):
            print(f'W A R N I N G :: {infiledir} not found !!!')
            continue
        if not os.path.isdir(bannerdir):
            print(f'W A R N I N G :: {bannerdir} not found !!!')
            continue

        jobFile = os.path.join(jobcardpath, handle+'.job') 
        shutil.copy(os.path.join(os.getcwd(),'sample.job.tmpl'), jobFile)
        os.chmod(jobFile, 0o755)
        mvaInFileName  = os.path.join(joboutpath, ptab_handle+'_mvaIn.root')
        mvaOutFileName = os.path.join(joboutpath, ptab_handle+'_mvaOut.root')
        histFile       = os.path.join(joboutpath, ptab_handle+'_'+args.tag+'_hist.root')
        runlogfile     = os.path.join(joboutpath, ptab_handle+'_'+args.tag+'.runlog')

        #bannerList = glob.glob(os.path.join(infiledir, '*.txt'))
        bannerList = [file for file in glob.glob(os.path.join(bannerdir, '*')) if '.txt' in file]
        if len(bannerList) == 0:
            print('No bannerFile found !!!')
            xsec = -1 if not isBkg else BkgDict.get(key)[0]
        else:
            banner = bannerList[0]
            xsecLine = [line for line in reversed(open(banner).readlines()) if '#  Integrated weight (pb)' in line]
            xsec = float(xsecLine[0].split(':')[-1].rstrip('\n')) 

        files    = [os.path.join(infiledir,item) for item in os.listdir(infiledir) if 'root' in item]
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
            #print(treeList)
            ttree = tfile.Get('Delphes')
            nEntries += ttree.GetEntries()

        print(f'\t\t cross-section : {xsec}, nEvents : {nEntries}')    
        ptab.add_row([ptab_handle, xsec, nEntries])

        BDT_WEIGHT1 = f'/home/gsaha/Work/DelphesML/Project2/BDT/IC_Model_ChiMass_{args.MX}/TMVAClassification_BDTD.weights.xml'

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
        elif args.tag == "BDT":
            replaceAll(jobFile, 'cutbased val', 'cutbased 0')
            replaceAll(jobFile, 'mvabased val', 'mvabased 1')
            replaceAll(jobFile, 'createMVATree val', 'createMVATree 1') if args.skim else replaceAll(jobFile, 'createMVATree val', 'createMVATree 0')
            replaceAll(jobFile, 'readMVA val', 'readMVA 1') if not args.skim else replaceAll(jobFile, 'readMVA val', 'readMVA 0') 
        else:
            replaceAll(jobFile, 'createMVATree val', 'createMVATree 0')
            replaceAll(jobFile, 'readMVA val', 'readMVA 0')            
            raise RuntimeError("Mention BDT or CUT")

        if MX == '60':
            replaceAll(jobFile, 'lowXmass val', 'lowXmass 1')
            replaceAll(jobFile, 'highXmass val', 'highXmass 0')
            replaceAll(jobFile, 'bdtThreshold score', 'bdtThreshold -0.03')
        elif MX == '90':
            replaceAll(jobFile, 'lowXmass val', 'lowXmass 0')
            replaceAll(jobFile, 'highXmass val', 'highXmass 1')
            replaceAll(jobFile, 'bdtThreshold score', 'bdtThreshold -0.01')

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

        if not isBkg:
            infodictSig[ptab_handle] = dict(sampleInfo())
            infodictSig[ptab_handle]['file']  = mvaInFileName
            infodictSig[ptab_handle]['xsec']  = xsec
            infodictSig[ptab_handle]['nEvents'] = nEntries
            sampledictSig.update(infodictSig)
        else:
            if key in ['Bkg_WZTo3LNu_012Jets', 'Bkg_ZZTo4L_012Jets']:
                infodictTrainBkg[ptab_handle] = dict(sampleInfo())
                infodictTrainBkg[ptab_handle]['file']  = mvaInFileName
                infodictTrainBkg[ptab_handle]['xsec']  = xsec
                infodictTrainBkg[ptab_handle]['nEvents'] = nEntries
                sampledictTrainBkg.update(infodictTrainBkg)
            else:
                infodictOtherBkg[ptab_handle] = dict(sampleInfo())
                infodictOtherBkg[ptab_handle]['file']  = mvaInFileName
                infodictOtherBkg[ptab_handle]['xsec']  = xsec
                infodictOtherBkg[ptab_handle]['nEvents'] = nEntries
                sampledictOtherBkg.update(infodictOtherBkg)
        
    procdict['SignalSamples'] = sampledictSig
    procdict['Train_Background'] = sampledictTrainBkg
    procdict['Other_Background'] = sampledictOtherBkg

    procdict['Lumi'] = 300000
    procdict['clsWeightFac'] = 1.0
    procdict['nEvPerSigSample'] = 20000
    procdict['outDir'] = '/home/gsaha/Work/DelphesML/DNN/Outputs'
    procdict['modelDir'] = '/home/gsaha/Work/DelphesML/DNN/Models'
    procdict['ftPlotDir'] = '/home/gsaha/Work/DelphesML/DNN/Outputs/FeaturePlots'
    

    with open(logfile, 'w') as logf:
        logf.write(str(ptab)+'\n')
    with open(cmdfile, 'w') as cmdf:
        for cmd in delphesCmdList:
            cmdf.write(cmd+'\n')

    #print(procdict)
    with open(yamlfile, 'w') as yamlf:
        yamldump = yaml.dump(procdict, yamlf, default_flow_style=False)

if __name__ == "__main__":
        main()
