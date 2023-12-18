list1 = [1,2,3]
list2 = ['a','b','c']

list3 = [[list1[position], list2[position]] for position in range(0,len(list1))]
print(list3)