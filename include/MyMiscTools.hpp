std::string fileExtension(std::string file){

    std::regex re(".*[^\\.]+\\.([^\\.]+$)");
    std::smatch result;
    if(std::regex_match(file,result,re))return result[1];
    else return "";

}

std::string fileNameWithoutExtension(std::string file){

    std::regex re("(.*[^\\.]+)\\.[^\\.]+$");
    std::smatch result;
    if(std::regex_match(file,result,re))return result[1];
    else return file;

}