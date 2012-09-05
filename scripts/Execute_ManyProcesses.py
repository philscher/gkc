import string
import sys
import time
import threading
import Queue
import os
import subprocess 
import copy

def getQueueFromList(list):
  queue = Queue.Queue()

  for p in list:
        task_unit = p
        queue.put(task_unit)
  return queue


result_list = []
#global list_results = []

class ProgramRun(threading.Thread):
    vlock = threading.Lock()
    def __init__(self, _startFunc, task_queue, result_list):
      threading.Thread.__init__(self)
      
      self.startFunc = _startFunc
      self.task_queue = task_queue 
      #self.result_list = copy.copy(result_list)

    def run(self):
      
      # Use infinite loop, not Queue.Full() because of threading
      while(True):

        try:
            # Get task unit, non-blocking
            task_unit = self.task_queue.get(False)
            print("Running task : " + str(task_unit) +  " Jobs  Left : " + str(self.task_queue.qsize()))
        except Queue.Empty:
            # Nothing to do anymore, return
            print("Nothing to do anymore .... exiting")
            return 

        ############################################################

        result = self.startFunc(task_unit)
       

        self.vlock.acquire()
        #self.result_list.append(result)
        result_list.append(result)
        self.vlock.release()

def ManyExecute(Func, par, nthreads):
    queue = getQueueFromList(par)
    #result_list = []
    threads = []
    for i in range(nthreads): 
        threads.append(ProgramRun(Func, queue, result_list))
        threads[-1].start()
    for i in range(len(threads)): threads[i].join()

    return result_list


