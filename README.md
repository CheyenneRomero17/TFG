# TFG: Análisis longitudinal de la atrofia cerebral hipocampal en estudios genéticos de la Enfermedad del Alzheimer

En este repositorio se pueden encontrar los distintos ficheros de código usados para realizar mi Trabajo de Fin de Grado titulado *Análisis longitudinal de la atrofia cerebral hipocampal en estudios genéticos de la Enfermedad del Alzheimer*.
La información médica usada en este proyecto está extraída del proyecto ADNI (Alzheimer’s Disease Neuroimaging Initiative, https://adni.loni.usc.edu/about/) y el *data management* de dicha información se ha realizado a partir de dos scripts de Python; *data management.ipynb* y *SNPs data.ipynb*. Por otro lado, la modelización de los datos y la visualización de los resultados de dichos modelos se han hecho a partir de scripts de R; *data modeling.R*, *model-to-xlsx.R* y *xlsx-to-fig.R*.

###### Python scripts
Así pues, los scripts usados para el *data management* de los datos de los participantes de este proyecto han sido
   - *data management.ipynb*: Para la visualización, limpieza y análisis de los datos generales de los participantes (sexo, edad, volumen hipocampal, etc) y
   - *SNPs data.ipynb*: Para la limpieza y el análisis de los datos genéticos de cada participante.
   
###### R scripts
Por otro lado, los scripts usados para la modelización y visualización de los modelos conseguidos han sido:
   - *data modeling.R*: Recoge los datos necesarios para realizar los modelos de regresión extrayéndolos de diversas tablas.
   - *model-to-xlsx.R*: Realiza los modelos de regresión lineales mixtos y guarda en una tabla los resultados.
   - *xlsx-to-fig.R*: Transforma las tablas de resultados a *heatmaps* para facilitar la comprensión de los resultados.

### Conclusiones
En este trabajo se han usado 12 variantes genéticas o SNPs de *risk* y 13 de *proxy* como variables explicativas para poder estudiar la atrofia cerebral en 13 subregiones del hipocampo de individuos a quienes se les ha realizado un seguimiento de la Enfermedad del Alzheimer desde su entrada en el proyecto ADNI hasta la actualidad. Para hacerlo se ha realizado un análisis longitudinal a través de modelos de regresión mixtos.
Los resultados obtenidos de los modelos muestran una influencia de la carga genética de la Enfermedad del Alzheimer en la atrofia de las subregiones del hipocampo a lo largo del tiempo. Esta influencia genética afecta de forma desigual a la atrofia hipocampal en función de como de avanzada se encuentre la enfermedad (diagnóstico clínico). No obstante, no se ha encontrado ninguna influencia conjunta que acelere el efecto de esta predisposición genética en la atrofia hipocampal.



Este trabajo de fin de grado ha sido supervisado por Natalia Vilor-Tejedor.
