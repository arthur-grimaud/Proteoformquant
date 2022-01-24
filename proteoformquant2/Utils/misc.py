import math
from itertools import tee

def truncate(num, n):
    integer = int(num * (10**n))/(10**n)
    return integer



# def match(mz1, mz2):
#     if abs( (1-mz1/mz2) * 1000000) <= 2:
#         return True
#     return False

# def match():
#     iter_2, last_match = tee(iter(l2))
#     d = {}
#     for i in l1:
#         d.setdefault(i, [])
#         found = False
#         while True:
#             j = next(iter_2, None)
#             if j is None:
#                 break
#             if matching(i, j):
#                 d[i].append(j)
#                 if not found:
#                     iter_2, last_match = tee(iter_2)
#                     found = True
#             else:
#                 if found:
#                     break
#         iter_2, last_match = tee(last_match)
#     return d