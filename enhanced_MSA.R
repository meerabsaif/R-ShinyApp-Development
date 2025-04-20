#Building an Enhanced version of Multiple Sequence Alignment App 
#Dividing into components for easy developing:
#1. Two methods to input sequences  (text & file)
#2. Alignment methods 
#3. Displaying results with consensus sequence
#4. Download option

library(shiny)                      #to build web application
library(msa)                        #provides tools for alignment & algorithms for msa
library(Biostrings)                 #handles biological sequences, including FASTA format
library(seqinr)                     #for reading and writing sequences
library(shinythemes)                #for improved UI themes

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Multiple Sequence Alignment App"),
  
sidebarLayout(                                            #divides the page into left sidebar and main panel 
  sidebarPanel(                                           #left panel, where user puts input
    width = 4,
  
#Selecting sequence type
radioButtons("seqType", "Choose Sequence Type:",               #radioButtons() creates bullet choice options 
choices = c("DNA" = "dna", "Protein" = "protein"),             #c("label" = "value")
selected = "dna"),                                             #setting default option to dna
      
#Creating Input method tabs
tabsetPanel(
tabPanel("Text Input", textAreaInput("seq", "Enter Sequences in FASTA Format:",         #textAreaInput() gives a text box
value = ">seq1\nATGCTAGCTAGCTACGATCG\n>seq2\nATGCTAGCTGGCTACGATCA",                     #\n used to separate line
height = "250px",
width = "150%")),
      
tabPanel("Upload File", fileInput("fastaFile", "Upload FASTA File",                     #creating another tab to input files 
accept = c(".fasta", ".txt"))                                                           #this only accepts files in fasta and text formats 
  )
),
      
#Creating a drop down menu to select an alignment method 
selectInput("methods", "Choose Alignment Method:",                         #selectInput() creates dropdown menu
choices = c("ClustalW", "ClustalOmega", "Muscle"),                         #these three alignment algorithms come with msa package 
selected = "ClustalW"),
      
#creating a button to Align sequences 
actionButton("alignBtn", "Align Sequences", class = "btn-primary"),        #class option improves overall look of the button by adding blue colour
      
#creating a button to download results
downloadButton("downloadBtn", "Download Results")

),
    
mainPanel(                                                                #output goes here 
    width = 8,                                                            
    h3("Sequence Information"),                                           #h3 used to add section headings 
    verbatimTextOutput("seqInfo"),                                        #displays sequence information                                                                  
      
    h3("Alignment Result"),
      verbatimTextOutput("alignmentResult"),                              #verbatimTextOutput (displays alignment in text format)
    
    h3("Consensus Sequence"),
    verbatimTextOutput("consensusSeq")                                    #verbatimTextOutput (displays consensus seq in text format)
    
    
    )
  )
)

server <- function(input, output) {                       #handles the logic behind what goes in the app 
                                                          #input takes value from UI, output displays what goes on the UI
#reactive value to store results
#reactive value: update automatically whenever new alignment is performed and triggers updates in apps 
  
values <- reactiveValues(                                 #reactiveValues() creates a reactive container, to place results later on
#creating 3 placeholders, setting them to NULL at start meaning they are empty 
sequences = NULL,                                         #to store input sequences 
alignment = NULL,                                         #to store alignment results 
consensusSeq = NULL                                       #to store consensus sequence

)
  
#Parse input sequences (from text or file)
getInputSequences <- reactive({                   #creates a reactive function to read input sequences from the user(either from text or a file)
                                                  #code will re-run whenever we enter a new text or upload a new file 
#Checking if file is uploaded
if (!is.null(input$fastaFile)) {                  #if the input file in not null then 
  temp_file <- input$fastaFile$datapath           #input$fastaFile$datapath gives the temporary path of the uploaded file to the server 
  
  } else {                                        #if no file is uploaded then we take whatever is typed in text area (input$seq)

    temp_file <- tempfile()                       #tempfile() creates a new temporary file path
    writeLines(input$seq, temp_file)              #writeLines() writes the given text into that temporary file
  }
  
#we write both file upload and text by writing them both into a file and reading from it 
#because these two functions: readDNAStringSet() and readAAStringSet() , read from a file and not directly from a string 
  
    
#Reading  sequences based on type
#these two functions come from BioStrings package
  
if (input$seqType == "dna") {                    #input$seqType comes from the radio button option in the UI ("dna" or "protein")
  sequences <- readDNAStringSet(temp_file)       #reads dna seq from fasta file
  } else {
    sequences <- readAAStringSet(temp_file)      #if seqtype isnt dna, then read aa seq 
 }

    return(sequences)                            #this returns parsed sequences 
                                                 #whenever we run getInputSequences(), we will get this value 
  })
  
#Displaying sequence information
output$seqInfo <- renderText({                         #output$seqInfo is linked to verbatimTextOutput("seqInfo")
  req(input$seq != "" || !is.null(input$fastaFile))    #req(), ensures this only runs if text or file is uploaded 
#text box is not empty [ input$seq != "" ]
#a file has been uploaded [ !is.null(input$fastaFile) ]
#if this condition isnt true, it wont run. This is to prevent errors
  
sequences <- getInputSequences()                        #getting sequences by calling getInputSequences() and storing them in the variable [ sequence ]
values$sequences <- sequences                           #saving it in the list of reactive values we created above, values.
  
  paste0(                                               #we are building a string here, paste0 () function combines different strings 
    "Number of sequences: ", length(sequences), "\n",
    "Sequence type: ", input$seqType, "\n",
    "Sequence lengths: \n",
    paste(
      names(sequences), ": ", width(sequences), " ",                    #gets names from fasta headers and length of each seq
      ifelse(input$seqType == "dna", "nucleotides", "amino acids"),     #ifelse() chooses either to write nucleotides or aminoa cids, based on seq type selection (dna or protein)
      collapse = "\n"                                                   #puts new line between each seqs info
    )
  )
})

#Performs alignment when button is clicked
observeEvent(input$alignBtn, {                         #observeEvent(), watches for when we click align sequences. 
req(input$seq != "" || !is.null(input$fastaFile))      #making sure we have input before running anything
  
sequences <- getInputSequences()                       #pulling user seqs again
values$sequences <- sequences                          #saving to salues$sequences, to reuse it
    
#performing alignment 
aligned <- msa(sequences, method = input$methods)      #msa() fn is run using selected method 
values$alignment <- aligned                            #storing aligned result in values
    
#Generating consensus sequence                         #consensus sequence are the most common residues at each aligned position
values$consensusSeq <- msaConsensusSequence(aligned)   #msaConsensusSequence, this generates consensus seq from the alignment, it is also stored to values 
    }, 
)

#Displaying alignment results
output$alignmentResult <- renderText({                 #output$alignmentResult is connected to verbatimTextOutput("alignmentResult") in the User interface
  req(values$alignment)                                #making sure alignment is computed 
  capture.output(print(values$alignment))              #print(values$alignment), prints alignment in a readable way
})                                                     #capture.output(), captures the printed text and returns as character string to display it

#Displaying consensus sequence                         #the comments written for alignment results, the logic applies exactly here 
output$consensusSeq <- renderText({
  req(values$consensusSeq)
  paste0("Consensus sequence: ", as.character(values$consensusSeq))     #as.character() fn converts consensus result to a plain string for displaying
  
})

#Downloading alignment results
output$downloadBtn <- downloadHandler(               #downloadHandler(), creates a download button
  filename = function() {                            #writing the name of the file that will be downloaded
    "alignment_result.txt"
},
  content = function(file) {                         #telling what to write in teh file when it is downloaded
    req(values$alignment)                            #to make sure alignment is available in the file before writing it
    
    writeLines(capture.output(print(values$alignment)), file)                   #writing alignment result to the file 
    writeLines(c("", "Consensus sequence:", as.character(values$consensusSeq)), 
        file, append = TRUE)                         #also inserting consensus seq in the file after the alignment result, adding a blank line between for layout
                                                     #append=true, to write new content after what is present in the file
  }                                                  #append=false, overwriting the file with new content
)}

shinyApp(ui, server)
