from matplotlib import pyplot as plt

## Função que retorna a sequência do arquivo fasta
def fastafile(file):

    fasta = open(file, 'r')
    fasta.readline()
    seq = ""
    for i in fasta.readlines():
        seq += i
    seq = seq.replace("\n","")
    fasta.close()
    return seq

def complement (string):
    # Função que retorna o complemento de uma sequência de DNA
    comple = ""
    for i in range (0, len(string)):
        if string[i] == "A":
            comple += "T"
        if string[i] == "C":
            comple += "G"
        if string[i] == "T":
            comple += "A"
        if string[i] == "G":
            comple += "C"
    return comple

'''def reading_frame(seq, number = 0):
    ini = seq.find("ATG", number)
    i=ini
    while (i+3<len(seq)):
        if seq[i:i+3] in ("TAG","TAA","TGA") :
            print(seq[ini:i+3])
            if (i%3 == ini %3):
                print(i)
                return i
        i+=3'''

def reading_frame(seq, number = 0):
    # Função que retorna o fim de um reading frame
    if number < len(seq):
        ini = seq.find("ATG", number)
        if ini > 0:
            i=ini
            while (i+3<len(seq)):
                if seq[i:i+3] in ("TAG","TAA","TGA") :
                    #print(seq[ini:i+3])
                    if (i%3 == ini %3):
                        #print(len(seq[number:i]))
                        return i
                i+=3
    return -1

def prom_box_value (seq1,seq2,box = 0):
    pont = 1
    match = 0
    mismatch = 0
    box10_value = 0
    box35_value = 0
    if box == 0 :
        for i in range (len(seq1)):
            if i==0:
                match = 0.8
            if i==1:
                match = 0.95
            if i==2:
                match = 0.45
            if i==3:
                match = 0.60
            if i==4:
                match = 0.50
            if i==5:
                match = 0.96
            mismatch = (1-match)/3
            if seq1[i]==seq2[i]:
                pont*=match
            else: pont*= mismatch
        box10_value = pont
        return box10_value ** (1/6)
    if box == 1:
        for i in range (len(seq1)):
            if i==0:
                match = 0.82
            if i==1:
                match = 0.84
            if i==2:
                match = 0.78
            if i==3:
                match = 0.65
            if i==4:
                match = 0.54
            if i==5:
                match = 0.45
            mismatch = (1-match)/3
            if seq1[i]==seq2[i]:
                pont*=match
            else: pont*= mismatch
        box35_value = pont
        return box35_value ** (1/6)

def prom_box (seq, seq_start):
    dummy = []
    for i in range (0,4):
        dummy.append(prom_box_value("TATAA", seq[seq_start-16+i:seq_start-10+i], 0))
    box10 = max(dummy)
    pos10_prom = seq_start-16+dummy.index(max(dummy))
    dummy.clear()
    for i in range (0,4):
        dummy.append(prom_box_value("TTGACA", seq[seq_start-40+i:seq_start-34+i], 1))
    box35 = max(dummy)
    pos35_prom = seq_start-40+dummy.index(max(dummy))
    return box10, pos10_prom, box35, pos35_prom


seq = fastafile("Mycoplasma genitalium.txt")
#seq = complement(seq)
#seq = seq[::-1]
#seq = complement(seq[::-1])

ini=[]
i = 0
ini.append((seq.find("ATG", 0)))
flag = 0

seqaux = seq + seq

box10 = []
pos10_prom = []
box35 = []
pos35_prom = []
a=0
b=0
c=0
d=0

while flag != -1:
    last = ini[-1] +1
    flag = seq.find("ATG", last)
    if flag != -1:
        ini.append(flag)


fim = [reading_frame(seq, i) for i in ini]
#fim=[]
#for i in ini:
#    if seq[i:i+3] in ("TAG","TAA","TGA"):
#        fim.append(i)

for index in  range (len(fim)):
    if fim[index] == -1:
        fim[index] = reading_frame(seqaux, ini[index]) - len(seq)



for value in ini:
    if value < 41:
        seqaux = seq[-41:] + seq
        a,b,c,d = prom_box(seqaux,value+41)
        box10.append(a)
        pos10_prom.append(b - len(seq))
        box35.append(c)
        pos35_prom.append(d- len(seq))
    else:
        a,b,c,d = prom_box(seq,value)
        box10.append(a)
        pos10_prom.append(b)
        box35.append(c)
        pos35_prom.append(d)


############################################################3333

seq = complement(seq)
ini2 = []
i = 0
ini2.append((seq.find("ATG", 0)))
flag = 0
seqaux = seq + seq

while flag != -1:
    last = ini2[-1] +1
    flag = seq.find("ATG", last)
    if flag != -1:
        ini2.append(flag)

fim2 = [reading_frame(seq, i) for i in ini2]
#fim2 = []
#for i in ini2:
#    if seq[i:i+3] in ("TAG","TAA","TGA"):
#        fim2.append(i)




for index in  range (len(fim2)):
    if fim2[index] == -1:
        fim2[index] = reading_frame(seqaux, ini2[index]) - len(seq)


box102 = []
pos10_prom2 = []
box352 = []
pos35_prom2 = []



for value in ini2:
    if value < 41:
        seqaux = seq[-41:] + seq
        a,b,c,d = prom_box(seqaux,value+41)
        box102.append(a)
        pos10_prom2.append(b - len(seq))
        box352.append(c)
        pos35_prom2.append(d- len(seq))
    else:
        a,b,c,d = prom_box(seq,value)
        box102.append(   a)
        pos10_prom2.append(b)
        box352.append(c)
        pos35_prom2.append(d)

'''
##############################################
'''
seq = seq[::-1]
i = 0
ini3=[]
ini3.append((seq.find("ATG", 0)))
flag = 0
seqaux = seq + seq

while flag != -1:
    last = ini3[-1] + 1
    flag = seq.find("ATG", last)
    if flag != -1:
        ini3.append(flag)

fim3 = [reading_frame(seq, i) for i in ini3]
#fim3 = []
#for i in ini3:
#    if seq[i:i+3] in ("TAG","TAA","TGA"):
#        fim3.append(i)


for index in range(len(fim3)):
    if fim3[index] == -1:
        fim3[index] = reading_frame(seqaux, ini3[index]) - len(seq)


box103 = []
pos10_prom3 = []
box353 = []
pos35_prom3 = []


for value in ini3:
    if value < 41:
        seqaux = seq[-41:] + seq
        a,b,c,d = prom_box(seqaux,value+41)
        box103.append(a)
        pos10_prom3.append(b - len(seq))
        box353.append(c)
        pos35_prom3.append(d- len(seq))
    else:
        a,b,c,d = prom_box(seq,value)
        box103.append(a)
        pos10_prom3.append(b)
        box353.append(c)
        pos35_prom3.append(d)



'''
####################################################
'''
seq = complement(seq[::-1])

i = 0
ini4=[]
ini4.append((seq.find("ATG", 0)))
flag = 0
seqaux = seq + seq

while flag != -1:
    last = ini4[-1] + 1
    flag = seq.find("ATG", last)
    if flag != -1:
        ini4.append(flag)

fim4 = [reading_frame(seq, i) for i in ini4]
#fim4 = []
#for i in ini4:
#    if seq[i:i+3] in ("TAG","TAA","TGA"):
#        fim4.append(i)


for index in range(len(fim4)):
    if fim4[index] == -1:
        fim4[index] = reading_frame(seqaux, ini4[index]) - len(seq)


box104 = []
pos10_prom4 = []
box354 = []
pos35_prom4 = []



for value in ini4:
    if value < 41:
        seqaux = seq[-41:] + seq
        a,b,c,d = prom_box(seqaux,value+41)
        box104.append(a)
        pos10_prom4.append(b - len(seq))
        box354.append(c)
        pos35_prom4.append(d- len(seq))
    else:
        a,b,c,d = prom_box(seq,value)
        box104.append(a)
        pos10_prom4.append(b)
        box354.append(c)
        pos35_prom4.append(d)

###################################################


ini.extend(ini2)
ini.extend(ini3)
ini.extend(ini4)
fim.extend(fim2)
fim.extend(fim3)
fim.extend(fim4)
box10.extend(box102)
box10.extend(box103)
box10.extend(box104)
box35.extend(box352)
box35.extend(box353)
box35.extend(box354)
pos10_prom.extend(pos10_prom2)
pos10_prom.extend(pos10_prom3)
pos10_prom.extend(pos10_prom4)
pos35_prom.extend(pos35_prom2)
pos35_prom.extend(pos35_prom3)
pos35_prom.extend(pos35_prom4)



cut = []
cut10 = []
cut35 = []
media10 = sum(box10)/len(box10)
media35 = sum(box35)/len(box35)
tamanho = []
indices = []
for x in range (len(ini)):
    cut.append(box35[x] * box10[x])
    if box35[x]*box10[x] > 0.1 and (fim[x]+3- ini[x]) > 30:
        indices.append(x)
        #cut.append(box35[x]*box10[x])
        cut10.append(box10[x])
        cut35.append(box35[x])
        i+=1
        tamanho.append(fim[i]+3-ini[i])




#print(indices)
#rint(box352)
#print(cut)




print("sequência de Mycoplasma genitalium: ")
print("      Start codon: ", ini[0])
print("      End codon: ", fim[0])
print("      score box 10:",box10[0])
print("      score box 35:",box35[0])
print("      score total:",box35[0]*box10[0])
print("      sequencia: ",seq[ini[0]:fim[0]+3])


#print(media)

plt.hist(box35, 50, color="darkslateblue")
plt.title('Histograma das pontuações das boxes 35')
plt.xlabel('pontuações das boxes 35')
plt.text(0.550, 55, f'Média = {media35:.4f}')
plt.show()
plt.cla()

plt.hist(cut35, 50, color="maroon")
plt.title('Histograma das pontuações das boxes 35')
plt.xlabel('pontuações das boxes 35')
plt.text(0.630, 111, f'Média = {sum(cut35)/len(cut35):.4f}')
plt.show()
plt.cla()

plt.hist(box10, 50, color="darkslateblue")
plt.title('Histograma das pontuações das boxes 10')
plt.xlabel('pontuações das boxes 10')
plt.text(0.6, 70, f'Média = {media10:.4f}')
plt.show()
plt.cla()

plt.hist(cut10, 50, color="maroon")
plt.title('Histograma das pontuações das boxes 10')
plt.xlabel('pontuações das boxes 10')
plt.text(0.425, 2.75, f'Média = {sum(cut10)/len(cut10):.4f}')
plt.show()
plt.cla()


plt.hist(cut, 50, color="goldenrod")
plt.title('Histograma das pontuações dos genes')
plt.xlabel('pontuações')
plt.text(0.30, 4000, f'Média = {sum(cut)/len(cut):.4f}')
plt.show()
plt.cla()


plt.hist(tamanho,50, color= "goldenrod")
plt.title('Histograma do tamanho dos genes')
plt.xlabel('Tamanho dos genes')
plt.text(2200, 700, f'Média = {sum(tamanho)/len(tamanho):.4f}')
plt.show()






'''
while i <  len(seq):
    try:
        a, b = reading_frame(seq, i)
        c,d,e,f = prom_box(seq,a)
    except TypeError:
        break
    ini.append(a)
    fim.append(b)
    box10.append(c)
    pos10_prom.append(d)
    box35.append(e)
    pos35_prom.append(f)
    i = a+1
'''
'''
if ("ATG" in seq[0:41]):
    seq1 = seq[len(seq)-41:len(seq)] + seq[0:41]
    find atg
    prom_box(seq1)
'''
'''
#for value in ini:
#    print(prom_box(ini,ini[value]))
print(len(ini))
#print(box10)

readingframe = {key:value for (key,value) in zip (ini,fim)}

#print(len(readingframe))
#print(ini)
#print(ini[-1])
#wtf = 0
#kekobreyer = []

while wtf != -1:
    last = ini[-1] +1
    wtf = seq.find("ATG", last)
    if wtf != -1:
        kekobreyer.append(wtf)
        ini.append(wtf)




print(kekobreyer)

seq1 = seq+seq
m = 0
n = 0

print(len(seq))
print(seq[len(seq)-1:len(seq)])

for x in kekobreyer:
    m,n = reading_frame(seq1,x)
    c,d,e,f = prom_box(seq,m)
    print(n)
    print(c,d,e,f)
    if n > len(seq):
        n = n - len(seq)
    fim.append(n)


print(len(fim))


'''



