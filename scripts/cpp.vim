"
"  Part of gkc project. 
"  Includes syntax highlitening for vim (www.vim.org)
"
"  for some defined types and other statements
"
"
"  Installation procedure :
"
"  copy to $HOME/.vim/after/syntax/
"
"


syn keyword cType cmplxd
syn keyword cType Array1z, Array2z, Array3z, Array4z, Array5z, Array6z
syn keyword	cStatement DMESG

" try to mark as red
syn keyword	cStatement omp_for
"syn keyword gkcSyntaxGroup transient
"hi gkcSyntaxGroup guifg=#ff00ff
"syn cluster GKCParaGroup add=gkcSyntaxGriyo

"syn keyword gkcSyntaxGroup omp_for





