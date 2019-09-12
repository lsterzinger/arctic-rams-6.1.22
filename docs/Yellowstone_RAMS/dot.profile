###My Personal Prompt
#PS1="\u@h:\w> "
PS1="\[\e[1;31;40m\]\u@\[\e[1;34;40m\]\H:\[e[1;32;40m\]\w>\[\e[0m\] "

export PS1

###Generic Aliases

alias rm='rm -i'
alias cp='cp -i'
alias mv='mv -i'
alias h="history | grep "
alias ls="ls --color"

### Other personal settings
export PATH=~/bin:.:$PATH
