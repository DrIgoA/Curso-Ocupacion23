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

- Conceptos básicos de estadística: Tipos de modelos y rol de los modelos en ciencia - Modelos estadísticos - Distribuciones de probabilidad - Probabilidad,
verosimilitud - Precisión, sesgo, exactitud. 
    - [Modulo 2 - Presentación](https://github.com/apgoijman/Curso-Ocupacion23/files/10883613/Modulo.2_Occupacion2023.pdf)
    - [Modulo 2 - Ejercicio MLE](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo2-Probabilidades/EjercicioModulo2.pdf)
    - [Modulo 2 - Script R para Ejercicio MLE](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo2-Probabilidades/Modulo2-EjemploLikelihhodBinomial.R)

- Introducción al enfoque Bayesiano: Comparación de inferencia frecuentista y Bayesiana - Inferencia Bayesiana y Teorema de Bayes -Componentes - Modelos jerárquicos - Ventajas y desventajas 
     - [Módulo 3 - Presentación](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/Modulo%203_Ocupacion2023.pdf) 

- Ejercicios con interfase R-JAGS para análisis con métodos Bayesianos (práctica) 
    - [Ejemplos Modulo 3 - Iniciación a R-JAGS](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo3-Bayes/Ejemplo%20Modulo3.md)
    - [Ejemplos Modulo 3 I](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo3-Bayes/Modulo%203%20-%20Ejemplo.R)
    - [Ejemplos Modulo 3 II](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo3-Bayes/Modulo3-%20Ejemplo%20II.R)
 

### DIA 2
- Modelos lineales generalizados con enfoque Bayesiano: Uso de interfase RJAGS para análisis Bayesianos – GLM – Introducción a efectos aleatorios y GLMM – GLMM Poisson y Binomial
    - [Módulo 4 - Presentación](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/Modulo%204_Occupacion2023.pdf)
- Ejercicios de GLM y GLMM con interfase R-JAGS para análisis con métodos Bayesianos (práctica) 
    - [Ejemplos Modulo 4, GLM y GLMM](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/Bayes/Modulo4-EjemploGLM-GLMM.R)
    - [Ejemplos Modulo 4, GLM Poisson](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/Bayes/Modulo4-EjemploPoissonGLM-Facu.R)

### DIA 3
- Detección imperfecta y modelos de ocupación para poblaciones: Solo una estación – Con covariables.
    - [Modulo 5 - Presentación](https://github.com/apgoijman/Curso-Ocupacion23/files/11214935/Modulo.5_Occupacion2023.pdf)
- Ejercicios de ocupación con interfase R-JAGS para análisis con métodos Bayesianos (práctica) y puesta en común.
- Trabajo práctico en proyectos finales (práctica)

[Formatos de bases de datos.pdf](https://github.com/apgoijman/Curso-Ocupacion23/files/10824363/Formatos.de.bases.de.datos.pdf)

### DIA 4
- Modelos de Ocupación de comunidades.

    - [Modulo 6 - Presentación](https://github.com/apgoijman/Curso-Ocupacion23/files/11254905/Modulo.6_Occupacion2023.pdf)
    - [Ejemplos Modulo 6 - Guía de un modelo de comunidades con DA](https://github.com/apgoijman/Curso-Ocupacion23/blob/main/EjerciciosEjemplos/Modulo6-OcupacionComunidad/Ejemplo%20Modulo%206.md)
    
- Ejercicios de ocupación con interfase R-JAGS para análisis con métodos Bayesianos y puesta en común (práctica).
- Trabajo práctico en proyectos finales (práctica)

### DIA 5
- Ejercicios de modelado, y trabajo en proyectos con datos propios.


## Modo de evaluación
La condición de aprobación será el 60% correspondiente a un 6 (seis). La evaluación del curso se realizará a través de un trabajo final que deberán entregar para su corrección 15 días después de finalizado el curso presencial. En la nota final se tendrá en cuenta la participación del estudiante a lo largo del curso y el trabajo final.

