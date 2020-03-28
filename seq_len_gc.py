#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

##导入模块，初始传递命令
import argparse
import re

parser = argparse.ArgumentParser(description = '该脚本用于统计 fasta 文件中，每条序列的长度和 GC 含量', add_help = False, usage = '\npython3.6 seq_len_gc.py -i [fasta_file] -o [output_file] -t [length or GC]\npython3.6 seq_len_gc.py --input [fasta_file] --output [output_file] --type [length or GC]')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-i', '--input', metavar = '[fasta_file]', help = '输入文件', required = True)
optional.add_argument('-o', '--output', metavar = '[output_file]', default = None, help = '输出文件，若缺省则直接打印在屏幕', required = False)
optional.add_argument('-t', '--type', metavar = '[length or GC]', default = 'length,GC', help = 'length，统计长度；GC，统计 GC 含量；length,GC，二者都统计（默认）', required = False)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

##读取并处理文件
#读取 fasta 文件并统计
seq_sum = {}
length_sum = 0; GC_sum = 0
with open(args.input, 'r') as seq_file:
	for line in seq_file:
		line = line.strip()
		if line[0] == '>':
			seq_id = line[1:len(line)]
			seq_sum[seq_id] = [0, 0]
		else:
			line_base = len(line)
			line_GC = len(re.findall('[GCgc]', line))
			seq_sum[seq_id][0] += line_base
			length_sum += line_base
			seq_sum[seq_id][1] += line_GC
			GC_sum += line_GC

seq_file.close()

#输出每条序列的长度或 GC 含量结果
if args.output:
	new_file = open(args.output, 'w')
	if 'length' in args.type and 'GC' in args.type:
		print('id\tlength\tGC', file = new_file)
		for key,value in seq_sum.items():
			print(f'{key}\t{value[0]}\t{round(100 * value[1] / value[0], 2)}', file = new_file)
	elif 'length' in args.type:
		print('id\tlength', file = new_file)
		for key,value in seq_sum.items():
			print(f'{key}\t{value[0]}', file = new_file)
	elif 'GC' in args.type:
		print('id\tGC', file = new_file)
		for key,value in seq_sum.items():
			print(f'{key}\t{round(100 * value[1] / value[0], 2)}', file = new_file)
	new_file.close()
else:
	if 'length' in args.type and 'GC' in args.type:
		print('id\tlength\tGC')
		for key,value in seq_sum.items():
			print(f'{key}\t{value[0]}\t{round(100 * value[1] / value[0], 2)}')
	elif 'length' in args.type:
		print('id\tlength')
		for key,value in seq_sum.items():
			print(f'{key}\t{value[0]}')
	elif 'GC' in args.type:
		print('id\tGC')
		for key,value in seq_sum.items():
			print(f'{key}\t{round(100 * value[1] / value[0], 2)}')

#总长总 GC（屏幕输出）
if 'length' in args.type and 'GC' in args.type:
	print(f'总碱基数： {length_sum}')
	print(f'GC 含量总百分比： {round(100 * GC_sum / length_sum, 2)}')
elif 'length' in args.type:
	print(f'总碱基数： {length_sum}')
elif 'GC' in args.type:
	print(f'GC 含量总百分比： {round(100 * GC_sum / length_sum, 2)}')
