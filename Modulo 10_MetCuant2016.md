

Módulo 10

MODELOS DE OCUPACIÓN DE

COMUNIDADES

Dra. Andrea P. Goijman

Curso de Posgrado: “Métodos cuantitativos de detección

imperfecta para el análisis de poblaciones y comunidades de

fauna silvestre”

Dpto. de Ciencias Naturales, UNRC

27 de Junio – 1 de Julio 2016





MODELOS DE COMUNIDADES

• Comunidad: ensamble de especies que

ocurren en un sitio

• Metacomunidad: Conjunto de

comunidades

• Riqueza de especies es la variable de

biodiversidad más ampliamente

utilizada.





MODELOS DE COMUNIDADES

Modelo jerarquico para comunidades o

meta comunidades

• Se puede modelar Riqueza de especies

• Modelos de ocupación de comunidades

(Dorazio & Royle 2005 - DR) - 3 niveles:

\1. Observación (error de medición)

\2. Especies (estado real)

\3. Comunidad





MODELOS DE COMUNIDADES

Modelo de comunidad como un

“hipermodelo” de distribución de un conjunto

de especies.

• Los parámetros para cada especie son

tratados como variables aleatorias – con

sus respectivas distribuciones *a priori*

• Los hiperparámetros de esas

distribuciones a priori, describen a la

comunidad.





MODELOS DE COMUNIDADES

Agregando complejidad se puede…

• Efectos de variables ambientales u

otras sobre la distribucion / abundancia

• Falsos negativos

• Dinamicas temporales

• Interacciones entre especies

• Niveles extra (ej. gremios)





MODELOS DE COMUNIDADES

• Detección imperfecta de especies puede

ser modelada con “Data augmentation”

– La inferencia sobre Riqueza total

inluye especies no registradas

• Inferencias a escala local (alfa), paisaje

(beta) y macroescala (gamma)





MODELOS DE COMUNIDADES

• Modelos que relacionan el numero de

especies observadas a covariables

(identidad de especies se pierde, no hay

error en Riqueza)

• N-mixtos (Detección imperfecta es igual

para todas las sp, no hay identidad de

especies)

• **Modelos de ocupación de**

**comunidades con efectos fijos de las**

**especies (sp. totalmente**

**independientes)**

• **Modelos de ocupación con las sp.**

**como efectos aleatorios (sólo**

**especies registradas).**





MODELOS DE COMUNIDADES

• Modelos de ocupación de comunidades con

efectos fijos de las especies (sp. totalmente

independientes)

– Lo mismo que ajustar un modelo para cada

sp por separado.

– Tiene la ventaja de poder comparar las

especies más fácilmente

– Especies poco frecuentes no mejoran su

estimación





MODELOS DE COMUNIDADES

• Modelos de ocupación con las sp. como

efectos aleatorios (sólo especies

registradas).

– Sp. siguen una distribución común de la

comunidad

– Mejora la precisión de las estimaciones para

especies poco frecuentes

– Para casos donde no nos interesan las

especies no registradas – la comunidad es

conocida

– Se puede modelar la correlación entre

ocupación y detectabilidad (abundancia

afecta detección)

– Covariables (heterogeneidad en deteccion





MODELOS DE COMUNIDADES

• Modelos de ocupación con las sp. como

efectos aleatorios (con DA).

– Lo mismo que en el modelo anterior,

incluyendo sp. que no fueron vistas, pero

sabemos que podrían ser parte de la

metacomunidad





REFERENCIAS

Marc Kery & J. Andy Royle. 2016. Applied

hierarchical modeling in ecology. Modeling

distribution, abundance and species

richness using R and BUGS. Volume 1:

Prelude and Static models. Academic

Press.

