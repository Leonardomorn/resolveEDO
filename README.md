Resolução de uma EDO - Equação diferencial ordinária - específica

Desenvolvedor:
Leonardo da Silva Camargo

Entrada e execução:
Para a execução do programa, basta chamar o Makefile e depois executar o arquivo executável resolveEDO. Ex: $make 
                                                                                                            $./resolveEDO

O programa resolve a EDO: 
     y´´+(x+1)y′−2y=(1−x²)e^(−x) , x∈(0,1), y(0)=−1 e y(1)=0

Primeiro método: Gera um sistema linear tri-diagonal com o método das diferenças finitas (MDF), sendo utilizado o método Gauss Seidel com vetores armazenar o Sistema Linear de forma a resolver o sistema linear.

Segundo método: Gera um sistema linear tri-diagonal com o método das diferenças finitas (MDF), sendo utilizado o método Gauss Seidel sem vetores armazenar o Sistema Linear de forma a resolver o sistema linear. Dessa forma, os vetores eram calculados antes.

Em cada um dos métodos é utilizado o número de passos do MDF como 5, 10 e 100.