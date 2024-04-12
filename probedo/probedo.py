# Author: Cameron F. Abrams <cfa22@drexel.edu>.
"""Probe PDB files

"""
import argparse as ap
import textwrap
import logging
import importlib.metadata
__probedo_version__ = importlib.metadata.version("probedo")
from .stringthings import banner, banner_message,aa_olc,join_ri,seqstr
from pidibble.pdbparse import PDBParser

logger=logging.getLogger(__name__)

class residue:
    def __init__(self,record):
        self.resName=record.resName
        self.seqNum=record.seqNum
        self.chainID=record.chainID
        self.iCode=record.iCode
    def __eq__(self,other):
        return self.resName==other.resName and self.seqNum==other.seqNum and self.chainID==other.chainID and self.iCode==other.iCode

def seq(args):
    pdbdesig=args.PDBfile
    if pdbdesig.endswith('.pdb'):
        pdbdesig=pdbdesig[:-4]
    p=PDBParser(PDBcode=pdbdesig).parse().parsed
    print(f'{args.PDBfile}: {len(p["ATOM"])} atoms')
    # resName: THR; chainID: G; seqNum: 90; iCode:
    rlist=[]
    for a in p["ATOM"]:
        r=residue(a.residue)
        if not r in rlist:
            rlist.append(r)
    print(f'{len(rlist)} unique residues')
    chains={}
    for r in rlist:
        if not r.chainID in chains:
            chains[r.chainID]=[]
        chains[r.chainID].append(r)
    for c in chains:
        a=""
        n=[]
        for r in chains[c]:
            a+=aa_olc[r.resName]
            n.append(join_ri(r.seqNum,r.iCode))
        print(f'Chain {c}:')
        print(seqstr(a,n))

def cli():
    commands={
        'seq':seq,
    }
    helps={
        'seq':'generate resid sequence map',
    }
    descs={
        'seq':'Use this command to generate a nice sequence map',
    }

    parser=ap.ArgumentParser(description=textwrap.dedent(banner_message),formatter_class=ap.RawDescriptionHelpFormatter)
    subparsers=parser.add_subparsers()
    subparsers.required=True
    command_parsers={}
    for k in commands:
        command_parsers[k]=subparsers.add_parser(k,description=descs[k],help=helps[k],formatter_class=ap.RawDescriptionHelpFormatter)
        command_parsers[k].add_argument('PDBfile',type=str,default=None,help='input pdb file (can fetch if not found)')
        command_parsers[k].set_defaults(func=commands[k])
        command_parsers[k].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
        command_parsers[k].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')
        command_parsers[k].add_argument('--diag',type=str,default='pestifer_diagnostics.log',help='diagnostic log file')
    
    args=parser.parse_args()
    args.func(args)