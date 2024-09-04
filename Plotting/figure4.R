library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(cowplot)
library(stats)
library(dbscan)

PCAfile = "~/Documents/Marklab/popfigure4/PCAmatrix.txt"
PCAheader = "~/Documents/Marklab/popfigure4/PCAmatrix.txt_header.txt"

HGSVC=c("HG00512_assem","HG00513_assem","HG00514_assem","HG00731_assem","HG00732_assem","HG03125_assem","NA12878_assem","NA19238_assem","NA19239_assem","NA24385_assem")

table = read.table(PCAheader, header = F)


table$V4 = -table$V4
table$V2[table$V1 %in% HGSVC] <- "HGSVC"
table$V2[grepl("GW", table$V1)] <- "CPC"
table= table[which(table$V1 != "CHM13_assem"),]


# Adjust shape values and labels
shape_values <- c("0" = 16, "1" = 17, "2" = 15)  # 16: circle, 15: square, 17: triangle
shape_labels <- c("1" = "male", "2" = "female")  # Rename shapes, omit "0" in legend
table$V2 <- factor(table$V2, levels = c("AFR", "EAS", "AMR", "SAS", "EUR", "gtex","CPC", "HGSVC","assem"))

# Convert V3 to a factor with specified labels
table$V3 <- factor(table$V3, levels = c("0","1", "2"))
legend_shapes <- setNames(c(15, 17), c("male", "female"))
# Color settings based on V2
color_values <- c("AFR" = "red", "EAS" = "yellow", "AMR" = "brown", "SAS" = "green", "EUR" = "blue", "assem" = "black", "CPC" = "orange","gtex" = "grey", "HGSVC" = "violet")

table$V2 <- factor(table$V2, levels = c("AFR", "EAS", "AMR", "SAS", "EUR", "gtex","CPC", "HGSVC","assem"))

# Order the table by 'V2'
table <- table %>% 
  arrange(V2)
 
minimal_grid_theme <- theme_minimal() +
                      theme(legend.position = "left",
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.border = element_blank(),
                            axis.line = element_line(colour = "black"))
                            
p1 <- ggplot(table, aes(x = V4, y = V5, color = V2, shape = V3)) +
  geom_point(aes(group = ifelse(V2 == "assem", 2, 1)), size = 2) +

  scale_color_manual(values = color_values, labels = c("1kg_AFR", "1kg_EAS", "1kg_AMR", "1kg_SAS", "1kg_EUR", "GTEx", "CPC", "HGSVC","HPRC")) +
  scale_shape_manual(values = shape_values, breaks = c("1", "2"), labels = c("male", "female")) +
  theme_minimal() +
  minimal_grid_theme + 
  theme(legend.position = "left") +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size=15)) +
  labs(x = "PC1", y = "PC2", color = "Populations", shape = "Gender")

p2 <- ggplot(table, aes(x = V5, y = V6, color = V2, shape = V3)) +
  geom_point(aes(group = ifelse(V2 == "assem", 2, 1)), size = 2) +

  scale_color_manual(values = color_values, labels = c("AFR", "EAS", "AMR", "SAS", "EUR", "GTEx", "HPRC", "HGSVC","HPRC")) +
  scale_shape_manual(values = shape_values, breaks = c("1", "2"), labels = c("male", "female")) +
  theme_minimal() +
  minimal_grid_theme + 
  theme(axis.text = element_text(size = 15), legend.text = element_text(size=15)) +
  guides(color = FALSE, shape = FALSE) +  # Remove legends
  labs(x = "PC2", y = "PC3", color = "Populations", shape = "Gender")
# Print the plot
p3 <- ggplot(table, aes(x = V6, y = V7, color = V2, shape = V3)) +
  geom_point(aes(group = ifelse(V2 == "assem", 2, 1)), size = 2) +

  scale_color_manual(values = color_values, labels = c("AFR", "EAS", "AMR", "SAS", "EUR", "GTEx", "HPRC", "HGSVC","HPRC")) +
  scale_shape_manual(values = shape_values, breaks = c("1", "2"), labels = c("male", "female")) +
  theme_minimal() +
  minimal_grid_theme + 
  theme(axis.text = element_text(size = 15), legend.text = element_text(size=15)) +
  guides(color = FALSE, shape = FALSE) +  # Remove legends
  labs(x = "PC3", y = "PC4", color = "Populations", shape = "Gender")
  
grid.arrange(p1, p2, p3, ncol = 3, widths = c(3.8,3,3))






AFR = (379, 412, 371, 385, 363, 360, 386, 382, 371, 354, 388, 352, 364, 378, 399, 391, 392, 379, 408, 369, 398, 380, 387, 351, 373, 403, 413, 394, 371, 369, 372, 362, 348, 355, 339, 381, 407, 380, 404, 395, 365, 379, 382, 384, 334, 401, 355, 380, 381, 428, 379, 378, 380, 393, 351, 360, 376, 371, 374, 367, 343, 364, 352, 376, 406, 407, 385, 396, 378, 418, 406, 386, 377, 418, 391, 392, 380, 358, 379, 397, 360, 414, 338, 354, 397, 361, 393, 397, 397, 368, 395, 396, 377, 360, 375, 359, 358, 381, 341, 384, 355, 388, 375, 354, 374, 376, 374, 393, 377, 347, 398, 406, 402, 408, 380, 441, 393, 367, 407, 395, 378, 404, 399, 355, 353, 335, 393, 361, 389, 371, 393, 365, 403, 379, 398, 350, 391, 348, 363, 376, 366, 406, 421, 349, 376, 379, 369, 340, 380, 398, 384, 421, 383, 371, 423, 370, 369, 366, 389, 329, 380, 353, 339, 353, 344, 424, 335, 356, 372, 378, 365, 357, 408, 410, 365, 368, 370, 408, 421, 361, 357, 393, 382, 426, 374, 378, 372, 380, 338, 411, 390, 378, 395, 394, 377, 373, 388, 389, 379, 363, 346, 386, 356, 390, 391, 381, 367, 323, 369, 374, 375, 349, 383, 393, 402, 382, 390, 394, 365, 397, 379, 354, 390, 342, 366, 389, 351, 409, 372, 360, 376, 382, 392, 388, 379, 359, 367, 375, 404, 401, 349, 528, 373, 374, 380, 426, 342, 412, 412, 360, 392, 364, 417, 385, 385, 381, 370, 376, 362, 409, 378, 391, 382, 401, 378, 387, 347, 369, 371, 366, 384, 382, 381, 374, 371, 398, 343, 388, 368, 381, 385, 348, 350, 388, 375, 377, 386, 366, 333, 431, 426, 341, 391, 407, 334, 384, 393, 368, 368, 437, 380, 361, 354, 390, 356, 377, 372, 345, 368, 406, 340, 377, 386, 416, 368, 392, 420, 365, 369, 324, 373, 390, 406, 385, 407, 357, 420, 380, 355, 356, 375, 361, 366, 352, 358, 380, 394, 403, 360, 427, 378, 397, 411, 382, 380, 369, 370, 392, 407, 357, 359, 383, 396, 362, 373, 362, 412, 387, 367, 388, 392, 366, 356, 383, 365, 388, 352, 410, 346, 375, 402, 405, 374, 387, 375, 355, 375, 382, 411, 373, 395, 363, 362, 379, 412, 363, 363, 400, 369, 383, 354, 378, 398, 418, 350, 438, 396, 364, 382, 385, 413, 380, 361, 345, 382, 385, 389, 367, 388, 387, 395, 354, 363, 370, 343, 376, 364, 347, 346, 325, 366, 395, 370, 344, 369, 361, 346, 366, 412, 396, 359, 371, 365, 353, 378, 344, 393, 375, 383, 413, 398, 387, 376, 391, 411, 373, 387, 390, 358, 413, 387, 379, 396, 389, 394, 375, 400, 379, 343, 393, 336, 369, 401, 404, 369, 400, 396, 361, 389, 370, 362, 383, 353, 351, 368, 392, 392, 386, 357, 386, 414, 432, 360, 362, 348, 395, 368, 345, 412, 390, 388, 361, 375, 357, 381, 390, 376, 396, 399, 340, 394, 349, 331, 381, 409, 386, 382, 352, 337, 423, 392, 404, 367, 389, 350, 374, 393, 385, 409, 352, 394, 353, 368, 373, 384, 353, 405, 371, 374, 343, 379, 398, 420, 370, 369, 363, 409, 375, 391, 372, 389, 363, 375, 342, 364, 411, 360, 372, 402, 379, 373, 350, 345, 382, 371, 409, 374, 404, 412, 406, 420, 423, 400, 422, 405, 420, 375, 381, 383, 390, 355, 405, 436, 408, 376, 360, 396, 399, 375, 372, 391, 368, 362, 400, 373, 364, 348, 413, 390, 332, 385, 378, 386, 399, 403, 335, 349, 383, 384, 379, 409, 382, 367, 377, 366, 383, 368, 376, 374, 393, 421, 388, 417, 366, 367, 369, 398, 391, 359, 360, 356, 369, 432, 396, 395, 367, 389, 333, 391, 371, 395, 341, 375, 351, 386, 407, 355, 413, 366, 383, 381, 351, 389, 388, 369, 376, 396, 354, 394, 390, 409, 378, 387, 400, 370, 365, 388, 387, 355, 366, 369, 380, 383, 374, 392, 372, 372, 344, 365, 369, 348, 360, 373, 335, 393, 351, 371, 391, 391, 349, 369, 401, 395, 426, 360, 398)
EAS = (357, 326, 324, 341, 281, 318, 328, 319, 321, 347, 355, 328, 357, 335, 344, 388, 362, 356, 353, 358, 358, 327, 348, 325, 327, 282, 355, 345, 311, 330, 388, 338, 370, 367, 361, 363, 352, 307, 370, 332, 317, 365, 304, 371, 353, 369, 322, 316, 400, 341, 373, 333, 329, 311, 328, 354, 366, 369, 335, 296, 364, 373, 349, 290, 332, 359, 311, 354, 329, 338, 328, 338, 341, 339, 366, 361, 347, 331, 347, 345, 351, 331, 344, 325, 364, 330, 304, 331, 348, 287, 339, 350, 299, 385, 372, 350, 324, 334, 326, 327, 325, 342, 333, 359, 319, 367, 355, 329, 371, 324, 370, 365, 314, 314, 333, 312, 370, 340, 376, 362, 344, 376, 347, 389, 341, 316, 361, 342, 341, 362, 310, 316, 391, 332, 353, 355, 329, 350, 315, 304, 357, 322, 355, 337, 349, 342, 343, 351, 347, 324, 331, 346, 377, 369, 335, 325, 352, 364, 358, 347, 319, 303, 327, 375, 305, 329, 352, 328, 341, 378, 333, 346, 337, 363, 384, 345, 376, 355, 322, 355, 337, 354, 350, 320, 319, 346, 328, 343, 350, 368, 345, 349, 355, 340, 334, 328, 356, 307, 319, 337, 354, 333, 339, 320, 322, 306, 348, 342, 355, 314, 340, 359, 309, 355, 399, 359, 328, 328, 324, 347, 373, 361, 365, 346, 329, 340, 344, 342, 328, 353, 363, 324, 324, 330, 341, 334, 364, 344, 324, 288, 333, 356, 345, 325, 323, 346, 353, 344, 347, 360, 361, 338, 343, 291, 328, 351, 344, 349, 399, 337, 332, 316, 318, 312, 338, 384, 328, 350, 339, 354, 358, 339, 355, 524, 330, 327, 313, 380, 339, 318, 374, 356, 308, 346, 367, 375, 353, 330, 348, 355, 354, 354, 346, 361, 313, 327, 345, 325, 343, 344, 348, 388, 317, 376, 348, 348, 344, 354, 343, 333, 345, 364, 345, 344, 354, 348, 352, 305, 366, 376, 329, 322, 364, 300, 341, 326, 330, 339, 320, 392, 333, 334, 364, 343, 381, 346, 345, 349, 360, 327, 336, 357, 384, 339, 367, 354, 361, 341, 328, 341, 333, 358, 324, 329, 329, 356, 326, 344, 330, 336, 292, 351, 339, 319, 347, 354, 334, 332, 358, 317, 337, 360, 325, 353, 334, 351, 325, 332, 322, 360, 331, 383, 333, 343, 375, 332, 368, 343, 354, 306, 348, 332, 328, 343, 364, 355, 377, 348, 330, 335, 335, 362, 345, 345, 361, 315, 353, 372, 364, 326, 351, 352, 331, 332, 323, 363, 351, 350, 328, 307, 302, 327, 344, 355, 310, 344, 359, 340, 375, 315, 354, 341, 500, 367, 315, 361, 342, 334, 353, 339, 340, 339, 324, 348, 327, 304, 337, 329, 350, 318, 328, 354, 379, 312, 366, 346, 325, 332, 368, 331, 347, 348, 322, 344, 319, 372, 341, 334, 280, 335, 386, 378, 336, 347, 351, 377, 352, 305, 322, 310, 337, 322, 378, 331, 327, 341, 341, 324, 338, 350, 309, 353, 333, 368, 341, 331, 315, 331, 316, 357, 353, 337, 368, 355, 384, 375, 343, 313, 341, 307, 346, 345, 335)
EUR = (412, 357, 362, 369, 286, 317, 313, 315, 365, 372, 365, 344, 354, 347, 322, 359, 302, 372, 314, 311, 313, 333, 343, 337, 337, 348, 367, 330, 331, 293, 374, 325, 340, 336, 380, 298, 312, 341, 384, 336, 338, 345, 305, 341, 332, 327, 328, 352, 316, 307, 363, 355, 307, 382, 328, 366, 366, 318, 321, 325, 317, 343, 320, 336, 333, 383, 346, 369, 318, 350, 384, 398, 381, 341, 353, 311, 307, 326, 400, 358, 368, 381, 308, 347, 377, 344, 311, 341, 332, 369, 348, 432, 348, 304, 342, 338, 334, 348, 363, 357, 365, 350, 357, 333, 348, 337, 347, 372, 343, 519, 545, 369, 333, 339, 378, 372, 383, 344, 349, 339, 345, 334, 341, 348, 326, 371, 341, 344, 351, 347, 318, 396, 340, 359, 334, 333, 325, 350, 335, 290, 318, 337, 356, 320, 376, 365, 320, 327, 356, 361, 302, 333, 372, 323, 303, 368, 326, 345, 320, 312, 318, 346, 319, 354, 343, 341, 400, 325, 329, 331, 305, 349, 343, 312, 352, 320, 319, 331, 374, 326, 342, 334, 355, 305, 292, 365, 330, 359, 356, 325, 372, 316, 353, 374, 319, 358, 364, 352, 301, 365, 376, 319, 347, 344, 381, 333, 300, 315, 359, 328, 319, 339, 391, 363, 374, 367, 354, 316, 327, 335, 353, 354, 362, 339, 348, 338, 359, 337, 349, 335, 368, 333, 341, 334, 339, 347, 331, 348, 333, 389, 364, 388, 353, 375, 365, 363, 329, 336, 385, 343, 328, 369, 345, 341, 368, 335, 356, 383, 297, 333, 343, 415, 354, 362, 360, 359, 327, 375, 334, 316, 321, 351, 310, 360, 306, 355, 301, 344, 356, 306, 381, 329, 356, 386, 352, 366, 333, 340, 333, 340, 374, 349, 339, 336, 344, 403, 319, 325, 337, 321, 364, 329, 331, 319, 328, 361, 353, 320, 349, 347, 371, 358, 316, 337, 367, 338, 325, 344, 331, 323, 366, 330, 334, 350, 376, 320, 346, 364, 344, 322, 355, 344, 326, 361, 350, 359, 364, 328, 306, 365, 328, 329, 321, 355, 340, 317, 340, 364, 366, 330, 350, 335, 303, 334, 322, 380, 360, 368, 363, 310, 338, 357, 344, 375, 382, 317, 354, 319, 305, 342, 362, 328, 365, 353, 300, 301, 340, 330, 343, 350, 371, 317, 324, 327, 333, 330, 320, 358, 332, 357, 338, 349, 347, 328, 336, 363, 336, 337, 344, 350, 378, 314, 335, 346, 323, 314, 333, 365, 352, 329, 321, 343, 333, 369, 352, 360, 317, 333, 321, 344, 729, 350, 352, 361, 316, 358, 338, 336, 312, 336, 348, 333, 383, 365, 335, 345, 471, 310, 333, 374, 302, 347, 331, 340, 360, 365, 313, 372, 342, 350, 321, 306, 349, 314, 321, 322, 313, 324, 341, 334, 333, 326, 354, 338, 331, 367, 309, 343, 340, 338, 337, 341, 356, 381, 396, 376, 377, 343, 320, 354, 340, 373, 313, 346, 320, 349, 311, 327, 285, 341, 323, 339, 349, 331, 308, 310, 363, 371, 338, 312, 327, 338, 315, 361, 319, 363, 330, 347, 357, 381, 334, 308, 338, 328, 329, 335, 286, 376, 366, 340, 308, 381, 371, 368, 348)
AMR = (340, 320, 330, 332, 329, 347, 324, 330, 353, 352, 345, 318, 305, 286, 370, 296, 314, 348, 337, 349, 342, 341, 324, 339, 356, 374, 320, 325, 332, 339, 351, 375, 328, 329, 334, 342, 394, 312, 321, 322, 308, 345, 342, 327, 355, 278, 330, 320, 333, 351, 308, 314, 334, 320, 301, 354, 331, 323, 298, 328, 345, 364, 328, 323, 355, 346, 333, 316, 345, 297, 335, 355, 350, 344, 343, 337, 333, 279, 349, 306, 323, 356, 320, 344, 324, 291, 338, 317, 299, 336, 371, 351, 317, 324, 334, 356, 333, 307, 328, 355, 323, 361, 380, 346, 335, 297, 344, 371, 324, 327, 352, 331, 303, 335, 317, 315, 376, 322, 336, 355, 361, 314, 347, 324, 329, 378, 362, 312, 312, 317, 360, 337, 342, 337, 336, 289, 340, 348, 354, 304, 370, 414, 353, 378, 347, 331, 342, 357, 307, 309, 324, 391, 312, 360, 311, 352, 311, 299, 324, 336, 331, 366, 364, 323, 351, 308, 328, 363, 331, 356, 322, 351, 321, 365, 339, 343, 314, 329, 344, 370, 331, 356, 353, 343, 328, 344, 319, 368, 329, 350, 342, 330, 290, 350, 357, 349, 321, 327, 370, 356, 319, 320, 296, 389, 312, 334, 345, 368, 355, 317, 331, 354, 334, 356, 351, 326, 346, 324, 538, 340, 367, 322, 343, 308, 331, 333, 347, 328, 338, 318, 396, 301, 326, 351, 343, 328, 318, 319, 352, 331, 311, 359, 318, 333, 356, 321, 344, 318, 325, 353, 358, 364, 364, 285, 309, 339, 337, 305, 342, 335, 366, 341, 302, 343, 353, 364, 339, 355, 354, 314, 333, 344, 366, 341, 353, 360, 322, 353, 344, 322, 302, 310, 332, 311, 318, 366, 351, 290, 329, 338, 336, 315, 327, 358, 378, 329, 321, 338, 328, 372, 310, 324, 297, 320, 322, 335, 351, 337, 339, 378, 357, 300, 329, 341, 336, 307, 306, 337, 326, 333, 326, 308, 351, 361, 338, 321, 299, 338, 350, 352, 355, 306, 326, 317, 313, 334, 317, 320, 321, 328, 338, 319, 350, 365, 328, 364, 309, 296, 329, 327, 355, 317, 336)
SAS = (358, 342, 313, 361, 338, 340, 326, 359, 368, 379, 358, 411, 355, 342, 356, 363, 352, 345, 343, 336, 371, 355, 361, 353, 348, 345, 345, 366, 347, 337, 338, 369, 349, 334, 345, 313, 361, 329, 361, 336, 359, 293, 338, 345, 362, 345, 351, 389, 381, 371, 342, 341, 375, 365, 350, 349, 358, 319, 315, 373, 367, 346, 324, 357, 387, 397, 331, 379, 338, 377, 347, 360, 366, 327, 340, 328, 348, 352, 376, 364, 338, 370, 364, 361, 321, 346, 360, 353, 327, 391, 319, 365, 319, 347, 368, 349, 360, 387, 388, 364, 337, 351, 362, 327, 345, 375, 363, 379, 341, 318, 340, 377, 356, 361, 369, 384, 354, 332, 370, 320, 351, 344, 332, 332, 350, 336, 351, 352, 364, 376, 368, 331, 376, 341, 356, 313, 349, 326, 357, 335, 282, 351, 352, 312, 358, 317, 331, 338, 360, 312, 365, 312, 364, 342, 376, 398, 365, 345, 376, 309, 367, 375, 339, 349, 379, 359, 359, 347, 350, 335, 344, 347, 357, 357, 374, 350, 388, 359, 311, 325, 310, 394, 362, 335, 353, 337, 376, 370, 355, 346, 407, 379, 345, 340, 306, 366, 375, 372, 327, 329, 380, 335, 371, 348, 376, 357, 330, 339, 352, 320, 352, 361, 328, 326, 340, 320, 310, 375, 347, 356, 359, 344, 342, 395, 386, 347, 333, 375, 348, 373, 373, 362, 332, 347, 358, 330, 338, 368, 357, 347, 336, 380, 349, 354, 374, 347, 317, 348, 369, 341, 350, 352, 303, 343, 355, 349, 332, 343, 341, 377, 357, 355, 330, 333, 397, 334, 336, 371, 368, 364, 317, 369, 367, 369, 354, 374, 344, 326, 352, 380, 363, 334, 333, 385, 366, 364, 364, 345, 384, 336, 362, 314, 352, 348, 359, 331, 356, 382, 373, 365, 363, 362, 365, 333, 359, 318, 334, 333, 333, 359, 368, 328, 309, 374, 348, 374, 333, 365, 325, 357, 351, 334, 353, 360, 355, 343, 365, 366, 365, 315, 369, 395, 357, 364, 377, 338, 375, 340, 330, 370, 376, 377, 334, 373, 330, 323, 358, 335, 356, 315, 383, 373, 316, 380, 360, 336, 314, 366, 355, 337, 336, 362, 351, 314, 357, 389, 385, 335, 323, 330, 341, 357, 364, 329, 328, 337, 326, 381, 342, 365, 354, 334, 358, 331, 315, 337, 329, 383, 347, 376, 367, 348, 336, 340, 320, 367, 391, 357, 361, 324, 356, 380, 316, 332, 348, 339, 365, 307, 344, 345, 337, 366, 327, 351, 355, 361, 415, 360, 339, 383, 330, 334, 368, 358, 337, 327, 359, 337, 362, 371, 344, 346, 354, 344, 337, 381, 313, 352, 363, 367, 365, 308, 341, 346, 385, 379, 335, 367, 321, 350, 329, 346, 377, 348, 343, 369, 354, 355, 364, 360, 345, 351, 325, 372, 347, 343, 357, 337, 339, 382, 314, 295, 307, 364, 339, 319, 373, 351, 348, 363, 376, 325, 366, 327, 355, 378, 335, 352, 387, 378, 305, 350, 341, 334, 324, 355, 321, 384, 344, 327, 371, 374, 336, 341, 352, 315, 349, 366, 351, 346, 328, 363, 346, 374, 364, 327)

# Create a list of vectors and their names
data_list <- list(AFR = AFR, EAS = EAS, EUR = EUR, AMR = AMR, SAS = SAS)

# Define a function to filter out outliers and return the filtered values and outliers
filter_outliers <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  filtered <- x[x >= mean_x - 5 * sd_x & x <= mean_x + 5 * sd_x]
  outliers <- x[!(x %in% filtered)]
  list(filtered = filtered, outliers = outliers)
}

# Apply the filter function to each vector in the list
filtered_data_list <- lapply(data_list, filter_outliers)

# Extract the filtered values and outliers into separate lists
filtered_values <- lapply(filtered_data_list, `[[`, "filtered")
outliers_values <- lapply(filtered_data_list, `[[`, "outliers")

# Combine the filtered values into a single data frame
data <- bind_rows(lapply(names(filtered_values), function(name) {
  data.frame(value = filtered_values[[name]], group = name)
}))

# Plotting
ggplot(data, aes(x = group, y = value)) +
  geom_boxplot() +
  labs(title = "", x = "", y = "Total aggregate copy number in each genome") +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.text = element_text(size = 15),     # Increase x-axis category label size
    axis.title = element_text(size = 15),
    panel.background = element_rect(fill = "white"),  # Set panel background to white
    panel.grid.major = element_line(color = "gray", size = 0.5)  # Add major grid lines
  )

mean(c(filtered_values$AMR , filtered_values$SAS , filtered_values$EUR , filtered_values$EAS))



thefile = "~/Documents/Marklab/popfigure4/dup_alltypes_fst.txt"
table = read.table(thefile, header = F, sep = "\t")
table  <- table [order(table$V4), ]

# Add row numbers as a new column
table <- table %>%
  mutate(Row = 1:length(rownames(table)))

# Extract label from V1 before the first "_"
table <- table %>%
  mutate(Label = sub("_.*", "", V1))


# Plotting the scatter plot with labels for points where V4 > 0.4
ggplot(table, aes(x = Row, y = V4)) +
  geom_point() +
  geom_text_repel(
    data = subset(table, V4 > 0.4), 
    aes(label = Label), 
    size = 3,
    box.padding = 0.3,        # Increase the padding around the text
    point.padding = 0.3,      # Increase the padding around the points
    max.overlaps = Inf,       # Allow for more overlaps if necessary
    nudge_x = -1000, nudge_y = 0,            # Nudge labels to the side
    direction = "y",          # Repel in the vertical direction
    hjust = 0                # Left-align the labels
  ) +
  labs(title = "", x = "", y = "Fst-like value") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),  # Center the title and increase its size
    axis.text.x = element_text(size = 24),  # Increase x-axis label size
    axis.text.y = element_text(size = 24)   # Increase y-axis label size
  )
  


thefile = "~/Documents/Marklab/popfigure4/all_analyze_noNA"
data = read.table(thefile, header = F, sep = "\t") 
data = data[which(data$V9 >= 20 & data$V8 != 'NA'),]
data <- data %>%
  mutate(Label = sub(";.*", "", V3))
data <- data %>%
  distinct(Label, .keep_all = TRUE)
# Plotting the scatter plot with labels for points where V4 > 0.4
# Define the ggplot scatter plot
newdata = data$V8
result <- dbscan(data.frame(newdata),eps = 0.15, minPts = 10)
data$cluster <- as.factor(3 - result$cluster)

stats <- data %>%
  group_by(cluster) %>%
  summarize(
    mean_V8 = mean(V8),
    sd_V8 = sd(V8),
    mean_V10 = mean(V10),
    sd_V10 = sd(V10)
  )

print(stats)
x_limits <- c(0, 50)
# Create the scatter plot
scatter_plot <- ggplot(data, aes(x = V8, y = V10, color = cluster)) +
  geom_point() +
  geom_text_repel(
    data = subset(data, (V10 > 3) | (V8 > 5 & V10 > 2) | (V8 > 10 & V10 > 0.8) | V8 > 15), 
    aes(label = Label),
    size = 3,
    box.padding = 0.2,  # Padding around the text box
    point.padding = 0.2  # Padding around the points
  ) +
  labs(title = "") +
  theme_minimal() + 
  labs(title = "", x = "", y = "", color = "cluster") +  # Set legend title to "cluster"
  theme(legend.text = element_text(size = 15),
    legend.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14),  # Center the title and increase its size
    axis.text.x = element_text(size = 24),  # Increase x-axis label size
    axis.text.y = element_text(size = 24),  # Increase y-axis label size
    plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"),  # Ensure consistent margins
    legend.position = c(0.5, 0.85)  # Move legend to the left
  ) +
  ylim(0, 3.5) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_color_manual(name="",values = c("3" = "black", "2" = "blue", "1" = "red"), labels = c("peak1", "peak2", "outliers"))

# Create the histogram
hist_plot <- ggplot(data, aes(x = V8)) +
  geom_histogram(bins = 300) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  labs(title = "", x = "", y = "Density") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis tick labels
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.y = element_text(size = 12),  # Ensure y-axis labels match scatter plot
    plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm")  # Ensure consistent margins
  )

# Combine the plots vertically
combined_plot <- plot_grid(hist_plot, scatter_plot, ncol = 1, rel_heights = c(1, 2), align = 'v')

# Display the combined plot
print(combined_plot)


library(dplyr)
library(dbscan)


thefile = "~/Documents/Marklab/popfigure4/all_tajimas"
tajimas_t = read.table(thefile, header = F, sep = '\t')


V1_vector <- as.character(tajimas_t$V1)
# Initialize an empty list to act as a dictionary
tajimas <- list()

# Process each value of V2_vector
for (i in seq_along(V1_vector)) {
  V1_list <- unlist(strsplit(V1_vector[i], ";"))
  V2_value <- abs(as.numeric(tajimas_t$V2[i]))
  for (key in V1_list) {
    if (!is.na(key)) {
      if (is.null(tajimas[[key]])) {
        tajimas[[key]] <- c()
      }
      tajimas[[key]] <- c(tajimas[[key]], V2_value)
    }
  }
}
  
mean_tajimas <- lapply(tajimas, mean, na.rm = TRUE)  



thefile = "~/Documents/Marklab/popfigure4/dup_alltypes_fst.txt"
table = read.table(thefile, header = F, sep = "\t")

thefile = "~/Documents/Marklab/popfigure4/all_analyze_noNA"
data = read.table(thefile, header = F, sep = "\t") 
data = data[which(data$V9 >= 20),]
data <- data %>%
  mutate(Label = sub(";.*", "", V3))
data <- data %>%
  distinct(Label, .keep_all = TRUE)
# Plotting the scatter plot with labels for points where V4 > 0.4
# Define the ggplot scatter plot
newdata = data$V8
result <- dbscan(data.frame(newdata),eps = 0.15, minPts = 10)
data$cluster <- as.factor(3 - result$cluster)



V2_vector <- as.character(table$V2)
# Initialize an empty list to act as a dictionary
dict <- list()

# Process each value of V2_vector
for (i in seq_along(V2_vector)) {
  V2_list <- unlist(strsplit(V2_vector[i], ","))
  V4_value <- as.numeric(table$V4[i])
  for (key in V2_list) {
    if (!is.na(key)) {
      if (is.null(dict[[key]])) {
        dict[[key]] <- c()
      }
      dict[[key]] <- c(dict[[key]], V4_value)
    }
  }
}
  
mean_dict <- lapply(dict, mean, na.rm = TRUE)  


V3_vector <- as.character(data$V3)

result_vector <- vector("numeric", length = length(V3_vector))

# Process each value of V3_vector
for (i in seq_along(V3_vector)) {
  x <- V3_vector[i]
  if (is.na(x) || x == "") {
    result_vector[i] <- NA
  } else {
    elements <- unlist(strsplit(x, ";"))
    values <- unlist(lapply(elements, function(ele) mean_dict[[ele]]))
    values <- values[!is.na(values)]
    if (length(values) == 0) {
      result_vector[i] <- NA
    } else {
      result_vector[i] <- mean(values)
    }
  }
}

data$result <- result_vector




filtered_data <- data %>% filter(!is.na(result))

recent_fst = filtered_data[which(filtered_data$cluster == 2),]$result
ancient_fst = filtered_data[which(filtered_data$cluster == 1),]$result

# Combine the data and compute ranks
combined_data <- c(recent_fst, ancient_fst)
ranks <- rank(combined_data)

# Split the ranks back into recent and ancient groups
recent_fst_ranks <- ranks[1:length(recent_fst)]
ancient_fst_ranks <- ranks[(length(recent_fst) + 1):length(ranks)]

# Calculate mean ranks
mean_rank_recent <- mean(recent_fst_ranks)
mean_rank_ancient <- mean(ancient_fst_ranks)

# Display the mean ranks for comparison
cat("Mean Rank - Recent_fst:", mean_rank_recent, "\n")
cat("Mean Rank - Ancient_fst:", mean_rank_ancient, "\n")

# Perform the Wilcoxon rank-sum test
wilcox_test <- wilcox.test(recent_fst, ancient_fst)

# Display the results of the Wilcoxon rank-sum test
print(wilcox_test)




result_vector2 <- vector("numeric", length = length(V3_vector))

# Process each value of V3_vector
for (i in seq_along(V3_vector)) {
  x <- V3_vector[i]
  if (is.na(x) || x == "") {
    result_vector2[i] <- NA
  } else {
    elements <- unlist(strsplit(x, ";"))
    values <- unlist(lapply(elements, function(ele) mean_tajimas[[ele]]))
    values <- values[!is.na(values)]
    if (length(values) == 0) {
      result_vector2[i] <- NA
    } else {
      result_vector2[i] <- mean(values)
    }
  }
}

data$tajimas <- result_vector2


filtered_data <- data %>% filter(!is.na(tajimas))

recent_tajimas = filtered_data[which(filtered_data$cluster == 2),]$tajimas
ancient_tajimas = filtered_data[which(filtered_data$cluster == 1),]$tajimas

# Combine the data and compute ranks
combined_data <- c(recent_tajimas, ancient_tajimas)
ranks <- rank(combined_data)

# Split the ranks back into recent and ancient groups
recent_tajimas_ranks <- ranks[1:length(recent_tajimas)]
ancient_tajimas_ranks <- ranks[(length(recent_tajimas) + 1):length(ranks)]

# Calculate mean ranks
mean_rank_recent <- mean(recent_tajimas_ranks)
mean_rank_ancient <- mean(ancient_tajimas_ranks)

# Display the mean ranks for comparison
cat("Mean Rank - Recent_fst:", mean_rank_recent, "\n")
cat("Mean Rank - Ancient_fst:", mean_rank_ancient, "\n")

# Perform the Wilcoxon rank-sum test
wilcox_test <- wilcox.test(recent_tajimas, ancient_tajimas)

# Display the results of the Wilcoxon rank-sum test
print(wilcox_test)









thefile = "~/Documents/Marklab/popfigure4/popanova.tsv"
data = read.table(thefile, header = F, sep = "\t")
data$V7 = data$V4/pmax(1, data$V2)
data$V8 = data$V6/pmax(1, data$V2)
data$V9 =  data$V8 - data$V7


# Sort V7 and V9 separately in ascending order
data_sorted_V7 <- data[order(data$V7), ]
data_sorted_V7$Index <- 1:nrow(data_sorted_V7)
data_sorted_V7$Variable <- 'Super-population'

data_sorted_V9 <- data[order(data$V9), ]
data_sorted_V9$Index <- 1:nrow(data_sorted_V9)
data_sorted_V9$Variable <- 'Subpopulation'

# Combine sorted data into a long format
data_long <- rbind(
  data.frame(Index = data_sorted_V7$Index, Value = data_sorted_V7$V7, Variable = data_sorted_V7$Variable),
  data.frame(Index = data_sorted_V9$Index, Value = data_sorted_V9$V9, Variable = data_sorted_V9$Variable)
)

# Plot data
ggplot(data_long, aes(x = Index, y = Value, color = Variable)) +
  geom_point() +
  ggtitle("") +
  theme_minimal() +
  labs(x = "", y = "Explained population variance of allele-specific copy numbers")


toplist = order(data$V9,decreasing = TRUE)[1:as.integer(length(rownames(data))/20)]



thefile = "~/Documents/Marklab/popfigure4/recombineestimate.txt"
data = read.table(thefile, header = F, sep = "\t")

data$V8 = data$V5 - data$V4

corr_cnv <- cor.test(data$V1, data$V7, method = "spearman")
corr_cnv

# ANOVA for V8 ~ V6
corr_length <- cor.test(data$V1, data$V8, method = "spearman")
corr_length

data <- data %>% mutate(loglen = log10(V8))

# Define breaks based on the range of log-transformed V9 values
min_loglen <- min(data$loglen, na.rm = TRUE)
max_loglen <- max(data$loglen, na.rm = TRUE)
log_breaks <- seq(floor(min_loglen), ceiling(max_loglen), by = 1)
original_labels <- 10^(log_breaks-3)  # Convert log breaks back to original scale for labeling
original_labels <- format(original_labels, scientific = FALSE)

# Create the scatter plot
# Create the scatter plot
ggplot(data, aes(x = V6, y = V1, color = loglen)) +
  geom_point() +
  scale_color_gradientn(
    colors = c("blue", "green", "yellow", "red"),
    breaks = log_breaks,  # Set breaks to the log-transformed values
    labels = original_labels   # Set formatted labels to avoid scientific notation
  ) +
  scale_x_continuous(
    breaks = seq(0, max(data$V7, na.rm = TRUE), by = 2),
    minor_breaks = seq(0, max(data$V7, na.rm = TRUE), by = 1)
  ) +
  labs(
    x = "Max mean population copy number of involving genes",
    y = "Multi-allelic LD",
    color = "Locus length\n(Kb)"
  ) +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size =24),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    panel.background = element_rect(fill = "white"), 
    legend.position = c(0.85, 0.85),  legend.justification = c(1, 1) 
  )







