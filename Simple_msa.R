#Builidng a simple msa app
#Dividing into components for easy developing:Components: 
#1. method for inputing sequences 
#2. method for performing alignment 
#3. download result option 


library(shiny)                       #to build web application
library(msa)                         #provides tools for aligning & contains algorithm for msa
library(Biostrings)                  #to handle biological sequences, handles fasta format too
library(seqinr)                      #handles sequence input and output, also provides tools foe seq manipulation and analysis

ui<- fluidPage (
  titlePanel ("Simple MSA App"), 
  sidebarLayout(
    sidebarPanel (                   #input here 
    
#creating input area for sequences 
textAreaInput("seq", "Enter the Sequences in Fasta Format:",                  #textAreaInput() gives a text box
value = ">seq1\nATGCTAGCTAGCTACGATCG\n>seq2\nATGCTAGCTGGCTACGATCA",                               #\n used to seperate line 
height = "250px" ,
width = "150%"),

    
#Creating a drop down menu to select an alignment method 
selectInput("methods", "Choose a Suitable Alignment Method:" ,                #selectInput() creates dropdown menu
choices = c( "ClustalW", "ClustalOmega", "Muscle"),
selected = "ClustalW"),

#Creaing a button for Alignment
actionButton ("alignbttn", "Align Sequences"),

#Creating a button for downloading the result obtained
downloadButton ("downbttn", "Download the Result")
      
      ),
    
    mainPanel (                      #output here 
     
      verbatimTextOutput("alignmentResult"),                                 #output area to show text output
                                               
       #verbatimTextOutput (displays alignment in text format)
 
    )
 
  )
)

server<- function (input, output) {  #main logic goes here

#reactive value to store alignment results:
#reactive value: update automatically whenever new alignment is performed
  
alignmentResult <- reactiveVal(NULL)                               
#reactiveVal  : It stores data, and when changes it updates other parts of the app too. Null: stating at 0 
#alignmentResult : created a container (which is empty rn) to place alignment results later on, when performed.
#here, it is used to Display the raw alignment results & Download full alignment data.


#perform alignment when button is clicked 
observeEvent (input$alignbttn, {                                #watches when we click "Align Sequences" button, then executes alignment process
 
 #parsing sequence from fasta format
 temp_file <- tempfile()                                   #creates a temporary file needed for a moment t store fasta input
 writeLines (input$seq, temp_file)                   #writing input to file. Now fasta data exists from where biostrings can find it
 sequences <- readDNAStringSet(temp_file)                   #now biostrings can easily read the file, and create a proper DNAstringset object, format needed for alignment algorithms
 
 #perform alignment 
 aligned <- msa(sequences, method = input$methods)          #performs msa with chosen method and returns an alignment object
 
    
 #storing results 
 alignmentResult(aligned)                                 #updating reactive values with results
  
  })

#displaying alignment results
output$alignmentResult <- renderText ({
  req(alignmentResult())                                  #making sure alignment exists before rendering
  capture.output(print(alignmentResult()))                #converts alignment object to text for displaying
  
  })

  
#downloading the alignment   
 output$downbttn <- downloadHandler(
   filename = function(){
     "alignment_result.txt"
   },
   content = function (file) {
   req(alignmentResult())
   writeLines(capture.output(print(alignmentResult())), file)
   }
  ) 
  
} 

shinyApp(ui,server)  
