pal_secreted <- c(
  "Actively secreted to blood" = "#B24F52",
  "Secreted in female reproductive system" = "#2A9D8F",
  "Unknown secreted location" = "#7C7C7C",
  "Secreted in male reproductive system" = "#A9D18E",
  "Secreted to extracellular matrix" = "#A6CEE3",
  "Secreted to digestive system" = "#264653",
  "Secreted in brain" = "#F4A5AE",
  "Secreted to other locations" = "#B89B74",
  "Digestive system" = "#264653",
  "Intracellular and Membrane" = "#457B9D",
  "Not part of the secretome annotation" = "lightgrey"
)

disease_class <- c(
  "Global" = "#E5E5E5",  # Add a color for Global
  "Healthy" = "#C9B28F",
  "Cardiovascular" = "#ED936B",
  "Metabolic" = "#E0C59A",
  "Cancer" = "#919FC7",
  "Neurologic" = "#7DC0A6",
  "Psychiatric" = "#7DC0A6",
  "Infection" = "#F9DA56",
  "Autoimmune" = "#DA8EC0"
)

pal_platform <- 
  c(
      "Significant in Olink" = "#6FC2EE",
      "Significant in SomaScan" = "#D82F88",
      "Significant in both" = "#312F7C",
      "Not significant in both" = "#E5E5E5"
    )


class_order <- c("Global", "Healthy", "Cardiovascular", "Metabolic", "Cancer", "Neurologic", "Autoimmune", "Infection")
disease_class_order <- c("Healthy", "Cardiovascular", "Metabolic", "Cancer", "Neurologic", "Autoimmune", "Infection")

