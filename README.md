# Modelado y estimación de ocupación para poblaciones y comunidades de especies bajo enfoque Bayesiano.
### CCT CONICET Mendoza 24 - 28 Abril 2023 
<img src="https://github.com/apgoijman/Curso-Ocupacion23/blob/main/varios/Imagen4.png" width=7% height=7%>


## Docentes
- Dra. Andrea Paula Goijman. INTA EEA La Consulta, Mendoza. 
- Dr. Facundo Contreras. CONICET, FCEyN Universidad Nacional de Rio Cuarto, Córdoba. 
- Dra. Vanesa Serafini. CONICET, FCEyN Universidad Nacional de Rio Cuarto, Córdoba. 

<img src="https://github.com/apgoijman/Curso-Ocupacion23/blob/main/varios/Imagen2.png" width=17% height=17%>    <img src="https://github.com/apgoijman/Curso-Ocupacion23/blob/main/varios/Imagen3.png" width=17% height=17%>   <img src="https://github.com/apgoijman/Curso-Ocupacion23/blob/main/varios/logo_giepco.png" width=17% height=17%>    <img src="https://github.com/apgoijman/Curso-Ocupacion23/blob/main/images.png" width=9% height=9%>

## Horario
Lunes a viernes de 8.30-12.15hs y de 13.45-19hs

## Objetivos del curso
El objetivo del curso es el de proveer bases teórico-prácticas de métodos cuantitativos bajo un enfoque Bayesiano para estimar parámetros para el análisis de poblaciones y comunidades de fauna silvestre con modelos jerárquicos. Se hará foco en los métodos de ocupación, que tienen en cuenta la detección imperfecta de flora y fauna silvestre.

## Expectativas
- Breve presentación de los alumnos
- ¿Cuáles son sus expectativas? 
- Nuestras expectativas...

![image](https://user-images.githubusercontent.com/124918841/222793849-89917531-59a4-4047-93ab-1148d1030d38.png)


## Requisitos
- R última versión
    - Windows: https://cran.r-project.org/bin/windows/base/
     -  Ubuntu: https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html#installing-r
     - Mac: https://cran.r-project.org/bin/macosx/
- Rstudio https://posit.co/download/rstudio-desktop/
- JAGS https://sourceforge.net/projects/mcmc-jags/

## Programa

### DIA 1
- Conceptos básicos para estudios de poblaciones y comunidades de fauna silvestre: Algunas definiciones - Poblaciones cerradas y abiertas – Comunidades – observaciones en ecología
    - [Modulo 1 - **Presentación**](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/Modulo%201_Occupacion2023.pdf)
    - [Modulo 1 - Ejercicio](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo1-Intro/Modulo%201.R)

- Conceptos básicos de estadística: Tipos de modelos y rol de los modelos en ciencia - Modelos estadísticos - Distribuciones de probabilidad - Probabilidad,
verosimilitud - Precisión, sesgo, exactitud. 
    - [Modulo 2 - **Presentación**](https://github.com/apgoijman/Curso-Ocupacion23/files/10883613/Modulo.2_Occupacion2023.pdf)
    - [*Ejercicio* probabilidades - Modulo 2](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo2-Probabilidades/Ejercicio1_Modulo2_Probabilidades.R)
    - [*Ejercicio* MLE - Modulo 2](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo2-Probabilidades/Ejercicio2-Modulo2.pdf)
    - [Script R para *Ejercicio* MLE - Modulo 2](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo2-Probabilidades/Ejercicio2-Modulo2-LikelihhodBinomial.R)

- Introducción al enfoque Bayesiano: Comparación de inferencia frecuentista y Bayesiana - Inferencia Bayesiana y Teorema de Bayes -Componentes - Modelos jerárquicos - Ventajas y desventajas 
     - [Módulo 3 - **Presentación**](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/Modulo%203_Ocupacion2023.pdf) 
- Ejercicios y ejemplos con interfase R-JAGS para análisis con métodos Bayesianos (práctica) 
    - [Ejemplo guia Modulo 3 - Iniciación a R-JAGS](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo3-Bayes/Ejemplo%20Modulo3.md)
    - [*Ejercicio* Modulo 3 - calculo de medias normal zorzales](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo3-Bayes/Ejemplo1-Modulo3-Media.R)
    - [*Ejercicio* Modulo 3 - poisson y jerarquico poisson ranas](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo3-Bayes/Ejemplo2-Modulo3-Jerarquico.R)
    - [*Ejercicio* Modulo 3 - Halcones peregrinos](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo3-Bayes/Ejercicio1%20-%20Modulo%203.R)
    - [*Ejercicio* Modulo 3 - Arboles](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo3-Bayes/Ejercicio2-Modulo%203.R)
 

### DIA 2
- Modelos lineales generalizados con enfoque Bayesiano: Uso de interfase RJAGS para análisis Bayesianos – GLM – Introducción a efectos aleatorios y GLMM – GLMM Poisson y Binomial
    - [Módulo 4 - **Presentación**](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/Modulo%204_Occupacion2023.pdf)
- Ejercicios de GLM y GLMM con interfase R-JAGS para análisis con métodos Bayesianos (práctica) 
    - [Ejemplos Modulo 4, GLM y GLMM - Asimetria fluctuante](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo4-GLMGLMM/Ejemplo1_Modulo4.R) y [Base de datos](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo4-GLMGLMM/ejemplo_AF.csv) 
    - [*Ejercicio* Modulo 4, GLM Poisson - insectos](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo4-GLMGLMM/Ejercicio2_Modulo4_PoissonGLM.R) y [**Base de datos - FALTA**]
    - [*Ejercicio* Modulo 4, GLMM Poisson - insectos](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo4-GLMGLMM/Ejercicio1_Modulo4_PoissonGLMM.R) y [Base de datos](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo4-GLMGLMM/Datos_Ejercicio1_Modulo4_PoissonGLMM.csv)

### DIA 3
- Detección imperfecta y modelos de ocupación para poblaciones: Solo una estación – Con covariables.
    - [Modulo 5 - **Presentación**](https://github.com/apgoijman/Curso-Ocupacion23/files/11214935/Modulo.5_Occupacion2023.pdf)
- Ejemplos y Ejercicios de ocupación con interfase R-JAGS para análisis con métodos Bayesianos (práctica) y puesta en común.
    - [Ejemplo - Ocupacion 1 especie con datos simulados](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/Ejemplo1_Modulo5_OcupacionSimple_datos-simulados.R)
    - [*Ejercicio* Ocupacion datos reales "bluebug"](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/Ejercicio1_Modulo5_OcupacionSimple_datos-reales.R) y [Base de datos](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/bluebug.txt)
    - [Ejemplo Ocupacion simple con datos en forma de VECTOR](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/Ejemplo-vector.R) y [Base de datos](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/datos_cm.csv)
   - Ejemplos de bases de datos. Como armarlas desde Excel con tablas dinámicas y desde R. [Excel para tablas dinamicas](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/CAPTURAS_TABLA.xlsx) + [Ejemplos de formatos de bases de datos - cómo acomodarlas en matrices de ocupación](https://github.com/apgoijman/Curso-Ocupacion23/files/10824363/Formatos.de.bases.de.datos.pdf) + [Archivo en R](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/Formatos%20Base%20de%20datos.R) + [Matriz de datos](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/datos_MATRIZ.csv) + [Datos longitudinal](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/datos_LONG_COMPLETOS.csv)
    - [*Ejercicio* Ocupacion datos reales Jilguero](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/Ejercicio-Modulo5-OcupacionSICLFA.R). En estos enlaces pueden encontrarse los archivos para correr el script: [Base de datos Jilguero](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/SICFLA.csv) + [Base de datos de covariables](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/Cova_sicfla.csv) + [Base de datos covariables anuales](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/Cova_sicflayear.csv) + En este enlace estan los datos corridos completos [Datos corridos completos](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo5-OcupacionSimple/SICFLA.rda)
- Trabajo práctico en proyectos finales (práctica)


### DIA 4
- Modelos de Ocupación de comunidades.
    - [Modulo 6 - **Presentación**](https://github.com/apgoijman/Curso-Ocupacion23/files/11254905/Modulo.6_Occupacion2023.pdf)
    - [*Ejemplos* Modulo 6 - Guía de un modelo de comunidades con DA](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/Ejemplo%20Modulo%206.md)
    
- Ejercicios de ocupación con interfase R-JAGS para análisis con métodos Bayesianos y puesta en común (práctica).
    - En este ejemplo vamos a ver modelos de ocupación de comunidades con DA (data augmentation) con y sin covariables. Los datos son tomados del libro de Kery&Royle 2016). También hay dos ejercicios dentro de este ejemplo. [Ejemplo MBH - Comunidad con DA](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/Modulo%206%20-%20Comunidad%20I%20-%20MHB%20Bird%20survey.R) + [Modelos ya corridos DA sin covariables](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/outDA.rda) + [Modelos ya corridos DA con covariables](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/outDAcov.rda) + [Modelos ya corridos DA con covariables2](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/outDAcov2.rda) 
    - [Ejercicio comunidad ratones 1](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/Modulo%206%20-%20Ejercicio%20Comunidad%20Roedores1.R) + [Datos roedores](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/Roedores.RData) + [Modelos ya corridos sp](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/roedsp.rda) + [Modelos ya corridos comunidad](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/roedcomu.rda) + [Modelos ya corridos N](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/roedN.rda)
    - [Ejercicio comunidad ratones 2](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/Ejercicio_ratones-Modulo6.R) + [Base de datos roedores](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/datos_roedores.csv)
    - [Ejercicio aves chaco](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/chaco_Aves.R) + [BDaves](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/species_sitiosAves.csv) + [BD covas](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/covas_sitiosAves.csv) + [numero reps por sitio](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/puntos_por_sitioAves.csv)
- Trabajo práctico en proyectos finales (práctica)

### DIA 5
- Ejercicios de modelado, y trabajo en proyectos con datos propios.


## Modo de evaluación
La condición de aprobación será el 60% correspondiente a un 6 (seis). La evaluación del curso se realizará a través de un trabajo final que deberán entregar para su corrección 15 días después de finalizado el curso presencial. En la nota final se tendrá en cuenta la participación del estudiante a lo largo del curso y el trabajo final.

