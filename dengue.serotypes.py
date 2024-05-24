
#### Uso RStudio para correr scripts de python. 
 #Para ellos escribir python en la terminal para comenzar a trabajar
python
 
# El módulo os es de la biblioteca estándar de Python

import os #Contiene funciones para realizar operaciones relacionadas con archivos, directorios, procesos y variables de entorno, entre otras cosas
 
os.getcwd() #para saber donde estoy trabajando
#>>>/home/rossa
 
os.makedirs(programacion_unlu) #para crear un directorio
 
os.chdir(programacion_unlu)# Cambiar al directorio destino

#Descargar las secuencias FASTA de la base de daton NCBI (base de datos de nucleotidos) y guardarlas en la carpeta creada "programacion_unlu". Esto lo hice moviendo los archivos de descaga sin la terminal

#Sigo usando el modulo os para listar los archivos de mi directorio 
os.listdir() 
#>>> "['dengue_genome1.fna', 'dengue_genome2.fna']"

# Cargar los genomas desde archivos FASTA en mi carpeta de trabajo
pip install biopython #Hacer esto fuera de python por q pip no lo reconoce


>>>Collecting biopython
  Obtaining dependency information for biopython from https://files.pythonhosted.org/packages/2a/38/3b971995c8bb2fad0b9809a61c8099fb1b2e579236e4c94fb0797825e171/biopython-1.83-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata
  Downloading biopython-1.83-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (13 kB)
Requirement already satisfied: numpy in /home/rossa/miniconda3/lib/python3.11/site-packages (from biopython) (1.26.1)
Downloading biopython-1.83-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (3.1 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 3.1/3.1 MB 2.7 MB/s eta 0:00:00
Installing collected packages: biopython
Successfully installed biopython-1.83"

from Bio import SeqIO #from Bio import SeqIO: Importa la clase SeqIO del paquete Bio, que es parte de Biopython. SeqIO proporciona funciones para leer y escribir archivos de secuencia en diferentes formatos.
#cargo genome1
genome1 = SeqIO.read("dengue_genome1.fna", "fasta") #read(): Es una función del módulo SeqIO que se utiliza para leer una única secuencia desde un archivo. Toma dos argumentos: el nombre del archivo a leer y el formato del archivo de entrada
print(genome1)#Para ver si la variable se genero correctamente

#>>>ID: NC_001477.1
#Name: NC_001477.1
#Description: NC_001477.1 Dengue virus 1, complete genome
#Number of features: 0
#Seq('AGTTGTTAGTCTACGTGGACCGACAAGAACAGTTTCGAATCGGAAGCTTGCTTA...TCT')

#Cargo genome2
genome2 = SeqIO.read("dengue_genome2.fna", "fasta") 
print(genome2)#Para ver si la variable se genero correctamente

#>>>ID: NC_001474.2
#Name: NC_001474.2
#Description: NC_001474.2 Dengue virus 2, complete genome
#Number of features: 0
#Seq('AGTTGTTAGTCTACGTGGACCGACAAAGACAGATTCTTTGAGGGAGCTAAGCTC...TCT')
 

# Crear una lista de secuencias
sequences = [genome1.seq, genome2.seq] #Esto crea una lista que contiene las secuencias de genome1 y genome2. La lista resultante, sequences, contendrá las secuencias de ambas variables genome1 y genome2

#Checkeo si genere esa variable
print(sequences)
#>>> [Seq('AGTTGTTAGTCTACGTGGACCGACAAGAACAGTTTCGAATCGGAAGCTTGCTTA...TCT'), Seq('AGTTGTTAGTCTACGTGGACCGACAAAGACAGATTCTTTGAGGGAGCTAAGCTC...TCT')]

# Escribir las secuencias en un archivo temporal
temp_file = "temp.fasta"
from Bio.SeqRecord import SeqRecord #from Bio.SeqRecord import SeqRecord: Importa la clase SeqRecord del paquete Bio. SeqRecord representa un registro de secuencia que contiene la secuencia biológica y metadatos asociados.
SeqIO.write([SeqRecord(seq, id=f"genome{i+1}") for i, seq in enumerate(sequences)], temp_file, "fasta")
>>>2 #En este caso, el valor 2 significa que dos secuencias se escribieron correctamente en el archivo temporal. Esto coincide con el hecho de que tienes dos secuencias en la lista sequences, ya que cada una de ellas se representa como un objeto SeqRecord que se escribió en el archivo temporal.

# Realizar el alineamiento utilizando MUSCLE
import subprocess
command = "/home/rossa/programacion_unlu/muscle.linux_intel64 -align temp.fasta -output out.fasta"
process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = process.communicate()


# Print the output and error
#stdout representa la salida estándar del proceso. Esto incluirá cualquier salida que el programa Muscle produzca durante su ejecución, como mensajes de progreso, estadísticas de alineación, etc.
#stderr representa la salida de error estándar del proceso. Esto incluirá cualquier mensaje de error que el programa Muscle pueda producir si algo sale mal durante su ejecución.

print(stdout.decode())
print(stderr.decode()) # El comando se ejecuta correctamente pero por algun motivo envia la salida standard a la salida de error 
#>>>muscle 5.1.linux64 [12f0e2]  8.3Gb RAM, 8 cores
#Built Jan 13 2022 23:17:13
#(C) Copyright 2004-2021 Robert C. Edgar.
#https://drive5.com

#Input: 2 seqs, avg length 10729, max 10735

#00:00 12Mb   CPU has 8 cores, running 8 threads
#00:00 71Mb    100.0% Calc posteriors
#01:18 71Mb    100.0% UPGMA5


# Leer el resultado del alineamiento
from Bio import AlignIO #from Bio import AlignIO: Importa la clase AlignIO del paquete Bio. AlignIO proporciona funciones para leer y escribir archivos de alineación de secuencias en diferentes formatos.

alignment = AlignIO.read("out.fasta", "fasta")

# Obtener la longitud del alineamiento
alignment_length = alignment.get_alignment_length()

# Inicializar una lista para almacenar los valores de similitud
similarity_values = []

# Calcular la similitud en cada posición del alineamiento (ojo que tiene que tener una distribucion el ciclo for)

for i in range(alignment_length):
    column = alignment[:, i]
    similarity = 1 if all(x == column[0] for x in column) else 0
    similarity_values.append(similarity)


# Generar el gráfico
import matplotlib.pyplot as plt #import matplotlib.pyplot as plt: Importa el módulo pyplot de la biblioteca Matplotlib, una biblioteca de visualización en Python. matplotlib.pyplot se usa comúnmente para crear gráficos y visualizaciones.

# Crear una lista de posiciones en el alineamiento para el eje x
positions = range(alignment_length)

# Crear el gráfico
plt.figure(figsize=(10, 6))  # Tamaño del gráfico
plt.plot(positions, similarity_values, marker='o', linestyle='-')  # Graficar los valores de similarity_values
plt.xlabel('Posición en el alineamiento')  # Etiqueta del eje x
plt.ylabel('Similaridad')  # Etiqueta del eje y
plt.title('Similaridad entre las secuencias en el alineamiento')  # Título del gráfico
plt.yticks([0, 1])  # Mostrar solo valores enteros en el eje y
plt.grid(True)  # Activar la cuadrícula
plt.savefig("plot.png") # Mostrar el gráfico
