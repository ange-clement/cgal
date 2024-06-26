#!/usr/bin/python

import os, sys, subprocess, datetime, time, signal, getopt
import numpy as np
import matplotlib.pyplot as plt

def main(argv):

  inputdir=""
  outputdir=""
  commit_hash=""
  try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:c:')
  except getopt.GetoptError:
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-i":
      inputdir = arg
    elif opt == "-o":
      outputdir = arg
    elif opt == "-c":
      commit_hash = arg

  print("Generating performance charts for inputdir =", inputdir)
  print("Outputdir =", outputdir)
  print("Commit hash =", commit_hash)

  exit_codes = {
   0 : "VALID_SOLID_OUTPUT",
   1 : "INPUT_IS_INVALID",
   2 : "OUTPUT_IS_INVALID",
   3 : "SIGSEGV",
   4 : "SIGABRT",
   5 : "SIGFPE",
   6 : "TIMEOUT"
  }

  current_run_data = {
   "VALID_SOLID_OUTPUT" : 0,
   "INPUT_IS_INVALID" : 0,
   "OUTPUT_IS_INVALID" : 0,
   "SIGSEGV" : 0,
   "SIGABRT" : 0,
   "SIGFPE" : 0,
   "TIMEOUT" : 0
  }

  filenames_per_codes = {}
  for key in current_run_data :
   filenames_per_codes[key] = []

  num_input = 0
  for filename in os.listdir(inputdir) :
    print("filename = ", filename)

    f = open(os.path.join(inputdir,filename))
    status = f.readline().strip();
    current_run_data[status] += 1
    filenames_per_codes[status].append(filename.rstrip('.log'))
    num_input = num_input+1

  # sort current_run_data by value
  current_run_data = {k: v for k, v in sorted(current_run_data.items(), key=lambda item: item[1], reverse=True)}

  # update chart data files
  date_now = datetime.datetime.now()
  date = str(date_now.year) +"-"+ str(date_now.month) +"-"+ str(date_now.day) +" "+ str(date_now.hour) +"h"+ str(date_now.minute) +"mn"

  for key_filename in current_run_data:
    f = open(os.path.join(outputdir+"/charts_data", key_filename+".txt"), "a+")
    f.write(str(current_run_data[key_filename]) + " " + commit_hash + " " + date + "\n")

  print("chart data updated")

  # update .pdf chart
  chart = plt.figure(figsize=(10, 7))
  colormap = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan","b","palegreen", "peachpuff"]
  plt.gca().set_prop_cycle('color', colormap)
  plt.style.use('tableau-colorblind10')
  for key_filename in current_run_data:

    f = open(os.path.join(outputdir+"/charts_data", key_filename+".txt"), "r")
    lines = f.readlines()
    x_number_values = []
    y_number_values = []
    i = 0
    for line in lines :
      if i < (len(lines) - 10) :
        i=i+1
        continue

      i=i+1
      words = line.strip().split()
      x_number_values.append(words[1]+"\n"+words[2]+"\n"+words[3])
      y_number_values.append(int(words[0]))
    plt.plot(x_number_values, y_number_values, marker='o', label=key_filename+": "+str(current_run_data[key_filename]))

  plt.xlabel("Version", fontsize=14)
  plt.ylabel("# of mesh", fontsize=14)
  plt.tick_params(axis="both", labelsize=9)
  plt.title("Benchmarking on " + str(num_input) + " meshes", fontsize=15)
  plt.legend(loc='lower left', bbox_to_anchor= (1.01, 0.58), ncol=1,
  borderaxespad=0, frameon=False)

  date_for_filename = str(date_now.year) +"-"+ str(date_now.month) +"-"+ str(date_now.day) +"-"+ str(date_now.hour) +"h"+ str(date_now.minute) +"mn"
  chart_filename = os.path.join(outputdir+"/charts","benchmarking_version_"+commit_hash+"-"+date_for_filename+".pdf")
  if os.path.isfile(chart_filename) :
     os.remove(chart_filename)
  chart.savefig(chart_filename, bbox_inches="tight")
  plt.close(chart)

  print("Robustness charts have been generated")

  # dump filenames per codes
  log_dirname = os.path.join(outputdir, "logs/"+commit_hash+"-"+date_for_filename)
  if not os.path.exists(log_dirname):
    os.mkdir(log_dirname)
  for key in filenames_per_codes :
   file = open(os.path.join(log_dirname, key+".txt"), "w+")
   for filename in filenames_per_codes[key] :
     file.write(filename + "\n")
   file.close()

  sys.exit()

if __name__ == "__main__":
  main(sys.argv[1:])
