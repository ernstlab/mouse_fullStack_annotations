#!/bin/csh -f
#  overlap_genome_context.sh.sh.cmd
#
#  UGE job for overlap_genome_context.sh.sh built Fri Feb  4 12:09:33 PST 2022
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/project/ernst/havu73/source/mm10_annotations/overlap/overlap_genome_context.sh.sh.joblog.$JOB_ID
#$ -o /u/project/ernst/havu73/source/mm10_annotations/overlap/overlap_genome_context.sh.sh.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/project/ernst/havu73/source/mm10_annotations/overlap/overlap_genome_context.sh.sh
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
#
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 2
# #
#  Name of application for log
#$ -v QQAPP=job
#  Email address to notify
#$ -M havu73@mail
#  Notify at beginning and end of job
#$ -m n
#  Job is not rerunable
#$ -r n
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "job serial"
  set qqidir    = /u/project/ernst/havu73/source/mm10_annotations/overlap
  set qqjob     = overlap_genome_context.sh.sh
  set qqodir    = /u/project/ernst/havu73/source/mm10_annotations/overlap
  cd     /u/project/ernst/havu73/source/mm10_annotations/overlap
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for overlap_genome_context.sh.sh built Fri Feb  4 12:09:33 PST 2022"
  echo ""
  echo "  overlap_genome_context.sh.sh directory:"
  echo "    "/u/project/ernst/havu73/source/mm10_annotations/overlap
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "overlap_genome_context.sh.sh started on:   "` hostname -s `
  echo "overlap_genome_context.sh.sh started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
#  module load java/jre-1.8.0_281
  module load intel
#
  echo overlap_genome_context.sh.sh "" \>\& overlap_genome_context.sh.sh.output.$JOB_ID
  echo ""
  /usr/bin/time /u/project/ernst/havu73/source/mm10_annotations/overlap/overlap_repeat_maskers_oarc.sh  >& /u/project/ernst/havu73/source/mm10_annotations/overlap/overlap_genome_context.sh.sh.output.$JOB_ID
#
  echo ""
  echo "overlap_genome_context.sh.sh finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/project/ernst/havu73/source/mm10_annotations/overlap/overlap_genome_context.sh.sh.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/job.log.serial
  if (`wc -l /u/project/ernst/havu73/source/mm10_annotations/overlap/overlap_genome_context.sh.sh.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
	head -50 /u/project/ernst/havu73/source/mm10_annotations/overlap/overlap_genome_context.sh.sh.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
	echo " "  >> /u/local/apps/queue.logs/job.log.serial
	tail -10 /u/project/ernst/havu73/source/mm10_annotations/overlap/overlap_genome_context.sh.sh.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  else
	cat /u/project/ernst/havu73/source/mm10_annotations/overlap/overlap_genome_context.sh.sh.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  endif
  exit (0)
