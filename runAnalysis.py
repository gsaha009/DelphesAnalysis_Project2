import sys
import os
import glob
import argparse
import shutil
import fileinput
from subprocess import Popen, PIPE
from joblib import Parallel, delayed


def runShellCmd(cmdList):
    proc = Popen(cmdList, stdout=PIPE, stderr=PIPE)
    while True:
        output = proc.stdout.readline()
        if proc.poll() is not None:
            break
        if output:
            print(output.strip().decode("utf-8"))
    rc = proc.poll()


def main():
    p = argparse.ArgumentParser(description='runDelphesAnalysis')

    p.add_argument('--Proc', type=str, action='store', required=True, help='XZ or HZ')
    p.add_argument('--MX', type=str, action='store', required=True, help='')
    p.add_argument('--ncore', type=int, action='store', required=False, default=2, help='')
    p.add_argument('--debug', action='store_true', required=False, help='')
    
    args = p.parse_args()
    
    tag = args.Proc+'_ChiMass_'+str(args.MX)
    infilePath = os.path.join(os.getcwd(),'JobCards')
    delphesCmdFilePath = os.path.join(infilePath, tag+'_delphesCmds.log')
    delphesCmdFile = open(delphesCmdFilePath)
    cmds = delphesCmdFile.readlines()
    cmdList = [(cmd.rstrip('\n')).split(' ') for cmd in cmds]
    #print(cmdList)
    popenCmdList = []
    for cmds in cmdList:
        exe = cmds[0]
        job = cmds[1]
        log = cmds[3]
        logf = open(log, 'wb')
        #logf.write('Log :----------> '+job+'\n')
        #logf.flush()
        if exe.startswith('X'):
            print('ignoring : {}'.format(job))
            continue
        popenCmdList.append([exe,job,logf])
        print(f'{exe}: {job}: {log}')

    ncores = args.ncore
    if not args.debug:
        delphes_run = Parallel(n_jobs=ncores)(delayed(runShellCmd)([delphesCMD[0],delphesCMD[1]]) for delphesCMD in popenCmdList)
        print(delphes_run)
    else:
        print('\n >>---------> Remove [--debug] to run analysis <---------<<')

if __name__ == '__main__':
    main()
