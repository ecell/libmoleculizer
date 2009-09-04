from urllib2 import build_opener, HTTPCookieProcessor, Request
from urllib import urlencode
from BeautifulSoup import BeautifulSoup, BeautifulStoneSoup
import re

#specify the root URL of the wiki
url_root = 'http://www.yeastpheromonemodel.org'

#specify the output bng file name
outputFile = open('model.bngl', 'w')


#Build opener
opener = build_opener(HTTPCookieProcessor)


#Authenticate
f=opener.open(Request(url_root + "/index.php?title=Special:Userlogin",urlencode(dict(wpName="Wiki parser",wpPassword="MM3vPVR",action="submitlogin"))))
f.close()


#######Parameters Section

category_url = url_root + "/index.php?title=Category:Parameters_-_Yeast_Pheromone_Response_Model"
outputFile.write('begin parameters\n')
param_number = 0

#since category pages can only contain 200 links, we have to loop through all
#parameter category pages
while (category_url != None):
    print "Opening category page: " + category_url

    #Retrieve Parameters category page
    f=opener.open(category_url)
    doc = f.read()
    f.close()

    #Parse the HTML from string
    soup = BeautifulSoup(doc)

    #find the link to the next 200 members of the category, if such a link
    #exists
    #This is done by first finding the text "previous 200", and then moving
    #forward to it's next sibling - if the sibling is a link entitled "next 200"
    #then that link will take us to the next category page.  Store URL for
    #future use, or set equal to None if there isn't a link.
    
    previous_200_text = soup.find(text=re.compile("previous 200"))
    if previous_200_text == None:
        next_200_a = None
    else:
        next_200_a = previous_200_text.findNextSibling('a')

    if next_200_a == None:
        category_url = None
    else:
        if next_200_a.string == "next 200":
            category_url = url_root + next_200_a['href']
            category_url = re.sub('&amp;', '&', category_url)
        else:
            category_url = None
        
    #locate all the parameter page urls within the current category page
    param_url_list = []
    for table in soup.findAll('table'):
        for a_tag in table.findAll('a', href=True):
                #get page name from url insteat of from a_tag title because some
                #page names have 'bad' characters that appear differently in url
            page_name = a_tag['href']           
            page_name = page_name[(page_name.index('wiki/') + 5):]
            param_url_list.append(url_root + "/index.php?title=" + page_name + "&action=raw")

    
    for param_url in param_url_list:
        #Retrieve Page
        print param_url
        f=opener.open(param_url)

        #Read the HTML contents into a string of text
        doc = f.read()
        f.close()

            #modeling tags are bracketed with < and > instead of < and >
            #replace these characters so beautifulsoup can recognize model tags
        doc = re.sub('>','>',doc)
        doc = re.sub('<','<',doc)

        #Parse the HTML from string
        soup = BeautifulSoup(doc)


        #find each parameter declared in the page (although with
        #current conventions there should only be one per page)
        #NOTE: Beautifulsoup makes tag names all lowercase!
        for param in soup.findAll('modelparameter'):

            param_number = param_number + 1
 
                #convert param (which is a sequence) to a string
            param_string = " ".join(["%s" % k for k in param])

                #replace the superscript tags with "^( )"
            param_string = re.sub('<sup>', '^(', param_string)
            param_string = re.sub('</sup>', ')', param_string)

                #remove '[[' and ']]'
            param_string = re.sub(r'\[','',param_string)
            param_string = re.sub(r'\]','',param_string)


                #remove extra spaces
            param_string = re.sub(' ', '', param_string)

                #write the parameter declaration to the file
            outputFile.write('%d' % param_number + '\t' + param_string + '\n')
    
outputFile.write('end parameters\n')



#####Species Section

category_url = url_root + "/index.php?title=Category:Species_-_Yeast_Pheromone_Response_Model"
molecule_types = ('\n\n\nbegin molecule types\n')
seed_species = ('\n\n\nbegin seed species\n')
species_number = 0

#since category pages can only contain 200 links, we have to loop through all
#species category pages
while (category_url != None):
    print "Opening category page: " + category_url

    #Retrieve Parameters category page
    f=opener.open(category_url)
    doc = f.read()
    f.close()

    #Parse the HTML from string
    soup = BeautifulSoup(doc)

    #find the link to the next 200 members of the category, if such a link
    #exists
    #This is done by first finding the text "previous 200", and then moving
    #forward to it's next sibling - if the sibling is a link entitled "next 200"
    #then that link will take us to the next category page.  Store URL for
    #future use, or set equal to None if there isn't a link.

    previous_200_text = soup.find(text=re.compile("previous 200"))
    if previous_200_text == None:
        next_200_a = None
    else:
        next_200_a = previous_200_text.findNextSibling('a')

    if next_200_a == None:
        category_url = None
    else:
        if next_200_a.string == "next 200":
            category_url = url_root + next_200_a['href']
            category_url = re.sub('&amp;', '&', category_url)
        else:
            category_url = None
            

    #locate all the species page urls within the current category page
    species_url_list = []
    for table in soup.findAll('table'):
        for a_tag in table.findAll('a', href=True):
                #get page name from url insteat of from a_tag title because some
                #page names have 'bad' characters that appear differently in url
            page_name = a_tag['href']           
            page_name = page_name[(page_name.index('wiki/') + 5):]
            species_url_list.append(url_root + "/index.php?title=" + page_name + "&action=raw")

    for species_url in species_url_list:
        #Retrieve Page
        print species_url
        f=opener.open(species_url)

        #Read the HTML contents into a string of text
        doc = f.read()
        f.close()

            #modeling tags are bracketed with < and > instead of < and >
            #replace these characters so beautifulsoup can recognize model tags
        doc = re.sub('>','>',doc)
        doc = re.sub('<','<',doc)

        #Parse the HTML from string
        soup = BeautifulSoup(doc)

        #find each molecule type declared in the page (although with
        #current conventions there should only be one per page)
        for molec in soup.findAll("modelmoleculetype"):
            species_number = species_number + 1

                #convert molec (which is a sequence) to a string
            molec_string = " ".join(["%s" % k for k in molec])

                #remove extra spaces
            molec_string = re.sub(', ', ',', molec_string)

                #remove '[[' and ']]'
            molec_string = re.sub(r'\[','',molec_string)
            molec_string = re.sub(r'\]','',molec_string)

                #add species_number and newline character
            molecule_types = molecule_types + '%d' % species_number + '\t' + molec_string + '\n'


        #find each seed species declared in the page (although with
        #current conventions there should only be one per page)
        for species in soup.findAll("modelseedspecies"):

                #convert species (which is a sequence) to a string
            species_string = " ".join(["%s" % k for k in species])

                #remove extra spaces
            species_string = re.sub(', ', ',', species_string)

                #remove '[[' and ']]'
            species_string = re.sub(r'\[','',species_string)
            species_string = re.sub(r'\]','',species_string)


                #add species_number and newline character
            seed_species = seed_species + '%d' % species_number + '\t' + species_string + '\n'

molecule_types = molecule_types + 'end molecule types\n'
seed_species = seed_species + 'end seed species\n'
outputFile.write(molecule_types)
outputFile.write(seed_species)




#######Reactions Section

category_url = url_root + "/index.php?title=Category:Reactions_-_Yeast_Pheromone_Response_Model"
outputFile.write('\n\n\nbegin reaction rules\n')
rule_number = 0

#since category pages can only contain 200 links, we have to loop through all
#reactions category pages
while (category_url != None):

    print "Opening category page: " + category_url

    #Retrieve Parameters category page
    f=opener.open(category_url)
    doc = f.read()
    f.close()

    #Parse the HTML from string
    soup = BeautifulSoup(doc)

    #find the link to the next 200 members of the category, if such a link
    #exists
    #This is done by first finding the text "previous 200", and then moving
    #forward to it's next sibling - if the sibling is a link entitled "next 200"
    #then that link will take us to the next category page.  Store URL for
    #future use, or set equal to None if there isn't a link.

    previous_200_text = soup.find(text=re.compile("previous 200"))
    if previous_200_text == None:
        next_200_a = None
    else:
        next_200_a = previous_200_text.findNextSibling('a')

    if next_200_a == None:
        category_url = None
    else:
        if next_200_a.string == "next 200":
            category_url = url_root + next_200_a['href']
            category_url = re.sub('&amp;', '&', category_url)
        else:
            category_url = None


    #locate all the reaction page urls within the current category page
    rxn_url_list = []
    for table in soup.findAll('table'):
        for a_tag in table.findAll('a', href=True):
                #get page name from url insteat of from a_tag title because some
                #page names have 'bad' characters that appear differently in url
            page_name = a_tag['href']
            page_name = page_name[(page_name.index('wiki/') + 5):]
            rxn_url_list.append(url_root + "/index.php?title=" + page_name + "&action=raw")


    for rxn_url in rxn_url_list:

        #Retrieve Page
        print rxn_url
        f=opener.open(rxn_url)

        #parse the page title out of the url
        page_title = rxn_url[(rxn_url.index('title=')+6):rxn_url.index('&action=')]
        outputFile.write('\n#' + page_title + '\n\n')

        #Read the HTML contents into a string of text
        doc = f.read()
        f.close()


            #modeling tags are bracketed with < and > instead of < and >
            #replace these characters so beautifulsoup can recognize model tags
        doc = re.sub('>','>',doc)
        doc = re.sub('<','<',doc)

        #Parse the HTML from string
        soup = BeautifulSoup(doc)

        #for each full reaction definition (rxn_full div), parse
        #out the elements of the reaction
        for rxn in soup.findAll("modelrxnfull"):

            #find each reaction equation in the full reaction
            for eqn in rxn.findAll("modelrxnrule"):
                rule_number = rule_number + 1

                #convert reaction equation (which is a sequence) to a string
                rxn_string = " ".join(["%s" % k for k in eqn])

                #remove any newlines within the reaction equation
                rxn_string = re.sub('\n[\ ]*', '', rxn_string)

                #add a '\' and a newline to the end of the reaction equation
                rxn_string = rxn_string + ' \\ \n\t'


                #add parameters to next line
                for rxn_param in rxn.findAll("modelrxnparam"):
                    rxn_string = rxn_string + '\t' +  ' '.join(["%s" % k for k in rxn_param])

                #remove all spaces
                rxn_string = re.sub(' ', '', rxn_string)
                #add spaces on either side of "+" (when "+" not preceded by "!")
                rxn_string = re.sub('(?<!!)\+', ' + ', rxn_string)
                #add space before "\"
                rxn_string = re.sub(r'\\', r' \\', rxn_string)

                #add space before reaction arrow
                rxn_string = re.sub('(?<!<)->', ' ->', rxn_string)
                rxn_string = re.sub('<->', ' <->', rxn_string)

                #add "\", newline, and tabs
                rxn_string = re.sub('->', '-> \\ \n\t\t', rxn_string)

                #remove '[[' and ']]'
                rxn_string = re.sub(r'\[','',rxn_string)
                rxn_string = re.sub(r'\]','',rxn_string)

                rxn_string = '%d' % rule_number + '\t' + rxn_string + '\n\n'
                rxn_string = re.sub('\t + ', '\t', rxn_string)

                outputFile.write(rxn_string)

outputFile.write('end reaction rules\n')


#add bng code to generate the network and create and SBML file
outputFile.write('\n\n\n\n\n')
outputFile.write('# Generation of the species and reactions\n')
outputFile.write('# with pheromone concentration set to zero.\n')
outputFile.write('generate_network({overwrite=>1});\n')
outputFile.write('# Write unequilibrated model to xml (model.xml)\n')
outputFile.write('writeSBML({prefix=>\"model\"});\n')
