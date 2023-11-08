# -*- 从定：utf-8 -*-
#!/usr/bin/env python3
'''
Created on Dec 3
@Author: Chen Tian
'''

import sys,os,logging,click
logging.basicConfig(filename=os.path.basename(__file__).replace('.py','.log'),
                    format='%(asctime)s: %(name)s: %(levelname)s: %(message)s',level=logging.DEBUG,filemode='w')
logging.info(f"The command line is:\n\tpython3 {' '.join(sys.argv)}")
#condon_table={}
def BARCODE(key_value_file): 
	"""Simple program to change barcode fastq file with 1 base mismatch"""
#f=open(key_value_file,'r')
	condon_table={}
	for i in key_value_file:
		i=i.strip().split()
		key=i[0]
		vaule=i[1]
		condon_table[key]=vaule
	return condon_table

@click.command()
@click.option('-1','--key_value_file',help='key value barcode file',type=click.File('r'),required=True)    ###此处设置了文件类型，可读或者可写
@click.option('-2','--barcode_raw',help='barcode raw',type=click.File('r'),required=True)
@click.option('-o','--output',help='barcode change file',type=click.File('w'),required=True)
def main(key_value_file,barcode_raw,output): 
#f1=open(barcode_raw)
#f2=open(barcode_change,'w')
	"""Simple program to change barcode fastq file with 1 base mismatch"""
	CONDON=BARCODE(key_value_file)
	for l in barcode_raw.readlines():
		l=l.strip()
		if l in CONDON:
			output.write(l.replace(l,CONDON[l]+'\n'))
		else:
			output.write(l+'\n')
	barcode_raw.close()
if __name__ == '__main__':                                                                     ###方便其他人调用此脚本
    main()
