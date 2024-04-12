# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" banner and resid/insertion splitting/joining routines
"""
import pandas as pd
import logging
logger=logging.getLogger(__name__)
import os
# _ANGSTROM_='Ångström'

import importlib.metadata

aa_olc={
    "ALA":"A",
    "ARG":"R",
    "ASN":"N",
    "ASP":"D",
    "CYS":"C",
    "GLN":"Q",
    "GLU":"E",
    "GLY":"G",
    "HIS":"H",
    "HSD":"H",
    "HSE":"H",
    "ILE":"I",
    "LEU":"L",
    "LYS":"K",
    "MET":"M",
    "PHE":"F",
    "PRO":"P",
    "SER":"S",
    "THR":"T",
    "TRP":"W",
    "TYR":"Y",
    "VAL":"V"
}
__probedo_version__ = importlib.metadata.version("probedo")

banner_message="""
    UNDER DEVELOPMENT
    {} v. {}
    https://probedo.readthedocs.io/en/latest/

    Cameron F. Abrams <cfa22@drexel.edu>

    Supported in part by Grants GM100472, AI154071, 
    and AI178833 from the NIH
    """.format(__package__.title(),__probedo_version__)

def banner(logf):
    my_logger(banner_message,logf,fill=' ',just='<')

def my_logger(msg,logf,width=67,fill='*',sep=', ',just='^',frame=''):
    """A fancy logger
    
    Parameters
    ----------
    msg: str, list
       the message to be logged, either as a single string or a list of strings
    
    logf: function
       writer; e.g., print, f.write, etc.

    width: int, optional
       linelength in bytes

    fill: str, optional
       single character used to fill blank spaces

    sep: str, optional
       single character used in join calls

    just: str, optional
       format character
    """
    fmt=r'{'+r':'+fill+just+f'{width}'+r'}'
    if frame:
        ffmt=r'{'+r':'+frame+just+f'{width}'+r'}'
    ll=' ' if just in ['^','>'] else ''
    rr=' ' if just in ['^','<'] else ''
    if frame:
        logf(ffmt.format(frame))
    if type(msg)==list:
        rr=' ' if ' ' not in sep else ''
        lnlst=[]
        for tok in msg:
            ll_sofar=sum([len(x) for x in lnlst])
            test_ll=ll_sofar+len(tok)+len(sep)*(len(lnlst)+1)
            if test_ll>width-2*(1+len(sep)):
                outstr=ll+sep.join(lnlst)+sep+rr
                logf(fmt.format(outstr))
                lnlst=[tok]
            else:
                lnlst.append(tok)
        outstr=ll+sep.join(lnlst)+' '
        logf(fmt.format(outstr))
    elif type(msg)==pd.DataFrame:
        for ln in msg.to_string().split('\n'):
            outstr=ll+ln+rr
            logf(fmt.format(outstr))
    else:
        lns=msg.split('\n')
        for ln in lns:
            outstr=ll+ln+rr
            logf(fmt.format(outstr))
    if frame:
        logf(ffmt.format(frame))
            

def split_ri(ri):
    """A simple utility function for splitting the integer resid and
    1-byte insertion code out of a string resid-insertion code
    concatenation
    
    Parameters
    ----------
    ri: the supposed resid-insertion concatenation
    
    Returns
    -------
    tuple(int, str): the integer resid and the 1-byte insertion code or '' if none
    """
    if ri[-1].isdigit(): # there is no insertion code
        r=int(ri)
        i=''
    else:
        r=int(ri[:-1])
        i=ri[-1]
    return r,i

def join_ri(resseqnum,insertion):
    if insertion=='':
        return str(resseqnum)
    return f'{resseqnum}{insertion}'

def ri_range(val,split_chars=['-','#']):
    the_split=[val]
    for c in split_chars:
        the_splits=[x.split(c) for x in the_split]
        the_split=[]
        for s in the_splits:
            the_split.extend(s)
    return [split_ri(x) for x in the_split]

def seqstr(s,n,gap=5,ll=60):
    print(n[0],n[-1])
    assert len(s)==len(n)
    nlines=len(s)//ll
    if len(s)%ll>0:
        nlines+=1
    ss=""
    for l in range(nlines):
        ss+='\n'
        left=l*ll
        right=min((l+1)*ll,len(s))
        nn=[split_ri(x)[0] for x in n[left:right]]
        # digits={}
        acc=True
        d=1
        dlabs=[]
        while acc:
            acc=any([x>10**(d-1) for x in nn])
            if acc:
                dlabs.append("".join([str((x%(10**d))//(10**(d-1))) for x in nn]))
                d+=1
        fst=True
        for ls in dlabs[::-1]:
            if fst:
                ts=ls.replace('0',' ')
            else:
                ts=ls
            fst=False
            tts=""
            for k,cc in enumerate(ts):
                if not k%gap:
                    tts+=cc
                else:
                    tts+=" "
            ss+=tts+'\n'
        ss+=s[left:right]+'\n'
    return ss
