# Source = https://stackoverflow.com/questions/8241099/executing-tasks-in-parallel-in-python

import threading

def task1(i):
    print(i)
    pass
def task2(i, j):
    print(i)
    print(j)
    pass
def task3():
    pass
def task4():
    pass
def task5():
    pass
def task6():
    pass

def dep1():
    t1 = threading.Thread(target=task1, args = (1,))
    t2 = threading.Thread(target=task2, args = (1, 2))
    t3 = threading.Thread(target=task3)

    t1.start()
    t2.start()
    t3.start()

    t1.join()
    t2.join()
    t3.join()

def  dep2():
    t4 = threading.Thread(target=task4)
    t5 = threading.Thread(target=task5)

    t4.start()
    t5.start()

    t4.join()
    t5.join()

def dep3():
    d1 = threading.Thread(target=dep1)
    d2 = threading.Thread(target=dep2)

    d1.start()
    d2.start()

    d1.join()
    d2.join()

d3 = threading.Thread(target=dep3)
d3.start()
d3.join()
