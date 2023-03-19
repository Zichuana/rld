'''
The following are useless
'''
import stronghold as sh

# x = int(input())
# y = int(input())
# sum = 0
# for i in range(x, y+1):
#     if i%2!=0:
#         sum+=i
# print(sum)


# f = open('rosalind_ini5.txt')
# flag = 0
# while True:
#     line = f.readline()
#     flag = flag + 1
#     if line:
#         if flag % 2 == 0:
#             print(line[:-1])
#     else:
#         break
# f.close()

# str = input()
# res = {}
# for word in str.split(' '):
#     if word in res:res[word] += 1
#     else:
#         res[word] = 1
# for key, value in res.items():
#    print(key, value)

with open(r"C:\Users\Zichuana\Downloads\rosalind_lcsm.txt", "r") as f:
    data = f.readlines()
    print(data)
dna_list = []
mid = ''
for i in data:
    if i[0] == '>':
        print(mid)
        dna_list.append(mid)
        mid = ''
    else:
        mid += i[:-1]
dna_list = dna_list[1:]
print('!!!!', len(dna_list[0]))
s = []
while 1:
    flag = input()
    if flag == '':
        break
    else:
        sub = input()
        print(sub)
        sub = sub.replace(' ', '')
        s.append(sub)
print(sh.Finding_a_Shared_Motif(s))
