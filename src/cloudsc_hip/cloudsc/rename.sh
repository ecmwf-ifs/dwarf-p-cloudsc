for file in *.cu; do 
    mv -- "$file" "${file%.cu}.cpp"
done
