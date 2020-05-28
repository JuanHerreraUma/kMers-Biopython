from Bio import SeqIO

for secuencia in SeqIO.parse('Llactis_Genome.txt', 'fasta'):

    lectura = secuencia.seq

    def cmers(lectura):

        kFreq = {'AA': 0,'AC': 0,'AG': 0,'AT': 0,'CA': 0,'CC': 0,'CG': 0,'CT': 0,'GA': 0,'GC': 0,'GG': 0,'GT': 0,'TA': 0,'TC': 0,'TG': 0,'TT': 0}

        for i in range(0, len(lectura) - 2 + 1):

            kmer = lectura[i:i + 2]

            if kmer in kFreq:
                kFreq[kmer] += (1 / len(lectura))

            else:
                print("No pertenece al diccionario de nucleótidos")

        return kFreq

    rf = cmers(lectura)

    print(secuencia.id, rf)
