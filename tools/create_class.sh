#!/bin/bash

echo ""
echo "===================================================================="
echo "This script creates the basic code for a new class file for genesis."
echo "===================================================================="
echo ""
echo "Please enter the module name that you want the class to be in:"
echo "(this is simply the containing folder for the code files)"
echo ""

cd ../lib/
select module in "Create new module" `ls -d */ | sed 's/\///g'` "Cancel"; do
    if [ "$module" == "Cancel" ]; then
        echo "Aborted. Nothing done."
        exit
    fi
    if [ "$module" == "Create new module" ]; then
        echo ""
        echo "--------------------------------------------------------------------"
        echo ""

        echo "Please enter a name for the new module:"
        echo ""
        read module
        module=${module,,}
        if [ "${module}" != `echo ${module} | sed 's/[^a-z]//g'` ]; then
            echo "Invalid module name. Name can only contain lowercase letters."
            echo "Aborted. Nothing done."
            exit
        fi
        if [ -d ${module} ]; then
            echo "Info: Module name already exists, it is not a new one."
        fi
        module="${module}/"
    fi
    if [ "${module}" != "" ]; then
        break
    fi
done

echo ""
echo "--------------------------------------------------------------------"
echo ""

# remove trailing / from name
module=${module%/}
echo "You selected module '$module'"

echo ""
echo "Please enter a name for the class, space separated for composed names:"
echo "Example: 'hello world' will create HelloWorld class."
echo ""
read name

echo ""
echo "--------------------------------------------------------------------"
echo ""

# remove surplus spaces
name=`echo ${name} | sed -r 's/^ *'// | sed -r 's/ *$//' | sed -r 's/  */ /g'`

if [ "${name}" != "`echo ${name} | sed 's/[^a-zA-Z0-9 ]//g'`" ]; then
    echo "Invalid class name. Name can only contain letters, digits and spaces."
    echo "Aborted. Nothing done."
    exit
fi

file_base=${name,,}
if [ "${module}" == "`echo ${name,,} | awk '{print $1;}'`" ]; then
    echo "The class name starts with the prefix of the module name."
    echo "Do you want to keep this in the file name, or shorten it for brevity?"
    echo "(the class name itself will be unchanged)"

    select yn in "Keep" "Shorten"; do
        if [ "${yn}" == "Shorten" ]; then
            file_base=`echo ${file_base} | sed -r 's/^\w*\ *//'`
        fi
        if [ "${yn}" != "" ]; then
            break
        fi
    done

    echo ""
    echo "--------------------------------------------------------------------"
    echo ""

fi

# make class name CamelCased and file name with underscores
class_name=`echo $name | sed -r 's/ ([a-z])/\U\1/g' | sed -r 's/^([a-z])/\U\1/g'`
file_name=`echo ${file_base,,} | sed -r 's/ /_/g'`

echo "About to create class '$class_name' in these files:"
echo "Header: ${module}/${file_name}.hpp"
echo "Code  : ${module}/${file_name}.cpp"
echo ""

if [ -f ${module}/${file_name}.hpp ] || [ -f ${module}/${file_name}.cpp ]; then
    echo "Files already exist. Aborting."
    exit
fi

echo "Continue?"

select yn in "Yes" "No"; do
    if [ "${yn}" == "No" ]; then
        echo "Aborted. Nothing done."
        exit
    fi
    if [ "${yn}" != "" ]; then
        break
    fi
done

echo ""
echo "--------------------------------------------------------------------"
echo ""

cap_module=`echo $module  | sed -r 's/([a-z])/\U\1/g'`
cap_file=`echo $file_name | sed -r 's/([a-z])/\U\1/g;s/_//g'`

mkdir -p ${module}

cat ../tools/class_template.hpp | sed "s/###CLASSNAME###/${class_name}/g" | \
    sed "s/###MODULNAME###/${module}/g;s/###FILENAME###/${file_name}/g"   | \
    sed "s/###CAPMODULNAME###/${cap_module}/g;s/###CAPFILENAME###/${cap_file}/g" > \
    ${module}/${file_name}.hpp

cat ../tools/class_template.cpp | sed "s/###CLASSNAME###/${class_name}/g" | \
    sed "s/###MODULNAME###/${module}/g;s/###FILENAME###/${file_name}/g"   | \
    sed "s/###CAPMODULNAME###/${cap_module}/g;s/###CAPFILENAME###/${cap_file}/g" > \
    ${module}/${file_name}.cpp

echo "Done."
