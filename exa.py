import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

class LinearProgrammingSolver:
    def __init__(self, root):
        self.root = root
        self.root.title("Solución de Programación Lineal")
        
        # Actualizar el número de variables y restricciones
        self.num_variables = 2  # Puede ajustar este valor según se necesite
        self.num_constraints = 3  # Se aumenta a 3 restricciones
        
        self.method = tk.StringVar(value="Método Simplex")
        self.optimization_type = tk.StringVar(value="Maximizar")
        
        self.create_widgets()
        
    def create_widgets(self):
        # Marco de tipo de optimización
        optimization_frame = ttk.LabelFrame(self.root, text="Tipo de Optimización")
        optimization_frame.pack(padx=10, pady=5, fill="x")
        
        optimization_types = ["Maximizar", "Minimizar"]
        for i, opt_type in enumerate(optimization_types):
            rb = ttk.Radiobutton(optimization_frame, text=opt_type, variable=self.optimization_type, value=opt_type)
            rb.grid(row=0, column=i, padx=5, pady=5)
        
        # Marco de función objetivo
        obj_frame = ttk.LabelFrame(self.root, text="Función Objetivo (Z)")
        obj_frame.pack(padx=10, pady=5, fill="x")
        
        self.obj_coeffs = []
        for i in range(self.num_variables):
            lbl = ttk.Label(obj_frame, text=f"Coeficiente de x{i+1}:")
            lbl.grid(row=0, column=2*i, padx=5, pady=5)
            entry = ttk.Entry(obj_frame, width=5)
            entry.grid(row=0, column=2*i+1, padx=5, pady=5)
            self.obj_coeffs.append(entry)
        
        # Marco de restricciones
        constraints_frame = ttk.LabelFrame(self.root, text="Restricciones")
        constraints_frame.pack(padx=10, pady=5, fill="x")
        
        self.constraints_coeffs = []
        self.constraints_signs = []
        self.constraints_rhs = []
        
        for i in range(self.num_constraints):
            row_entries = []
            for j in range(self.num_variables):
                lbl = ttk.Label(constraints_frame, text=f"Coeficiente de x{j+1}:")
                lbl.grid(row=i, column=3*j, padx=5, pady=5)
                entry = ttk.Entry(constraints_frame, width=5)
                entry.grid(row=i, column=3*j+1, padx=5, pady=5)
                row_entries.append(entry)
            # Signo
            sign_cb = ttk.Combobox(constraints_frame, values=["<=", ">=", "="], width=3)
            sign_cb.grid(row=i, column=3*self.num_variables, padx=5, pady=5)
            sign_cb.current(0)
            self.constraints_signs.append(sign_cb)
            # Lado derecho
            rhs_entry = ttk.Entry(constraints_frame, width=5)
            rhs_entry.grid(row=i, column=3*self.num_variables+1, padx=5, pady=5)
            self.constraints_rhs.append(rhs_entry)
            self.constraints_coeffs.append(row_entries)
        
        # Marco de restricciones de variables
        var_constraints_frame = ttk.LabelFrame(self.root, text="Restricciones de Variables")
        var_constraints_frame.pack(padx=10, pady=5, fill="x")
        
        self.var_constraints = []
        for i in range(self.num_variables):
            lbl = ttk.Label(var_constraints_frame, text=f"x{i+1} ≥ 0")
            lbl.grid(row=0, column=i*2, padx=5, pady=5)
            var = tk.IntVar(value=1)
            cb = ttk.Checkbutton(var_constraints_frame, variable=var)
            cb.grid(row=0, column=i*2+1, padx=5, pady=5)
            self.var_constraints.append(var)
        
        # Marco de selección de método
        method_frame = ttk.LabelFrame(self.root, text="Método de Solución")
        method_frame.pack(padx=10, pady=5, fill="x")
        
        methods = ["Método Simplex", "Método de la M grande", "Método de las dos fases"]
        for i, method in enumerate(methods):
            rb = ttk.Radiobutton(method_frame, text=method, variable=self.method, value=method)
            rb.grid(row=0, column=i, padx=5, pady=5)
        
        # Botón de solución
        solve_button = ttk.Button(self.root, text="Solucionar", command=self.solve)
        solve_button.pack(pady=10)
        
        # Marco de salida
        output_frame = ttk.LabelFrame(self.root, text="Solución")
        output_frame.pack(padx=10, pady=5, fill="both", expand=True)
        
        self.output_text = tk.Text(output_frame, wrap="word")
        self.output_text.pack(fill="both", expand=True)
        
    def solve(self):
        # Recolectar datos
        try:
            c = [float(entry.get()) for entry in self.obj_coeffs]
            if self.optimization_type.get() == "Minimizar":
                c = [-coeff for coeff in c]
            A = []
            b = []
            signs = []
            for i in range(self.num_constraints):
                row = [float(entry.get()) for entry in self.constraints_coeffs[i]]
                A.append(row)
                b.append(float(self.constraints_rhs[i].get()))
                signs.append(self.constraints_signs[i].get())
            # Restricciones de variables
            var_constraints = [var.get() for var in self.var_constraints]
        except ValueError:
            messagebox.showerror("Error de Entrada", "Por favor, ingrese valores numéricos válidos.")
            return
        
        # Implementar el método seleccionado
        if self.method.get() == "Método Simplex":
            self.simplex_method(c, A, b, signs, var_constraints)
        elif self.method.get() == "Método de la M grande":
            self.big_m_method(c, A, b, signs, var_constraints)
        elif self.method.get() == "Método de las dos fases":
            self.two_phase_method(c, A, b, signs, var_constraints)
        else:
            messagebox.showerror("Error de Método", "Por favor, seleccione un método válido.")
            return
        
    def simplex_method(self, c, A, b, signs, var_constraints):
        self.output_text.delete(1.0, tk.END)
        
        # Verificar si hay restricciones '>=' o '='
        if '>=' in signs or '=' in signs:
            self.output_text.insert(tk.END, "El Método Simplex estándar no puede resolver restricciones '≥' o '='.\n")
            return
        
        # Convertir a forma estándar
        num_variables = len(c)
        num_constraints = len(A)
        
        # Extender coeficientes de la función objetivo
        c_extended = c[:]
        
        # Agregar variables de holgura
        tableau = []
        basis = []
        slack_var_index = num_variables
        for i in range(num_constraints):
            row = A[i][:]
            slack = [0] * num_constraints
            if signs[i] == "<=":
                slack[i] = 1
                c_extended.append(0)
                basis.append(slack_var_index)
                slack_var_index +=1
            else:
                # No se pueden manejar '>=' o '=' en el Método Simplex estándar
                self.output_text.insert(tk.END, "El Método Simplex estándar no puede resolver restricciones '≥' o '='.\n")
                return
            row.extend(slack)
            row.append(b[i])
            tableau.append(row)
        
        # Mostrar el tableau inicial
        self.output_text.insert(tk.END, "Tableau Inicial:\n")
        self.display_tableau(tableau, c_extended, basis)
        
        # Iniciar iteraciones
        iteration = 0
        max_iterations = 100  # Para evitar ciclos infinitos
        while iteration < max_iterations:
            iteration +=1
            self.output_text.insert(tk.END, f"\nIteración {iteration}:\n")
            # Calcular las ganancias relativas (cj - zj)
            zj = [0]*len(c_extended)
            for i in range(len(basis)):
                for j in range(len(c_extended)):
                    zj[j] += c_extended[basis[i]] * tableau[i][j]
            cj_zj = [c_extended[j] - zj[j] for j in range(len(c_extended))]
            
            self.output_text.insert(tk.END, f"Ganancias relativas (cj - zj): {['{0:.2f}'.format(num) for num in cj_zj]}\n")
            
            # Verificar optimalidad
            if all(value <= 1e-5 for value in cj_zj[:-1]):
                self.output_text.insert(tk.END, "Se encontró la solución óptima.\n")
                break
            
            # Variable entrante (mayor cj - zj)
            entering_candidates = [j for j in range(len(cj_zj)-1) if cj_zj[j] > 1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes. El problema puede tener soluciones múltiples.\n")
                break
            entering = max(entering_candidates, key=lambda j: cj_zj[j])
            self.output_text.insert(tk.END, f"Variable entrante: x{entering+1}\n")
            
            # Calcular razones (b_i / a_i)
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))
            self.output_text.insert(tk.END, f"Razones: {['{0:.2f}'.format(num) if num != float('inf') else 'Inf' for num in ratios]}\n")
            
            # Verificar si es ilimitado
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "La solución es ilimitada.\n")
                return
            
            # Variable saliente (mínima razón positiva)
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: x{basis[leaving]+1}\n")
            
            # Pivoteo
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            # Actualizar base
            basis[leaving] = entering
            
            # Mostrar tableau
            self.display_tableau(tableau, c_extended, basis)
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones. El problema puede no tener solución óptima.\n")
            return
        
        # Extraer solución
        solution = [0]*len(c_extended)
        for i in range(len(basis)):
            if basis[i] < len(solution):
                solution[basis[i]] = tableau[i][-1]
        
        z = sum(c_extended[i]*solution[i] for i in range(len(c_extended)))
        if self.optimization_type.get() == "Minimizar":
            z = -z
        self.output_text.insert(tk.END, f"\nSolución Óptima:\n")
        for i in range(num_variables):
            self.output_text.insert(tk.END, f"x{i+1} = {solution[i]:.2f}\n")
        self.output_text.insert(tk.END, f"Z = {z:.2f}\n")
        
    def display_tableau(self, tableau, c_extended, basis):
        # Mostrar el tableau de forma formateada
        header = "Base\t"
        for j in range(len(c_extended)):
            header += f"x{j+1}\t"
        header += "LD\n"
        self.output_text.insert(tk.END, header)
        for i in range(len(tableau)):
            row_str = f"x{basis[i]+1}\t"
            for value in tableau[i]:
                row_str += f"{value:.2f}\t"
            self.output_text.insert(tk.END, row_str + "\n")
        
    def big_m_method(self, c, A, b, signs, var_constraints):
        self.output_text.delete(1.0, tk.END)
        
        M = 1e6  # Valor grande para M
        num_variables = len(c)
        num_constraints = len(A)
        
        c_extended = c[:]
        tableau = []
        basis = []
        artificial_vars = []
        var_index = num_variables
        
        for i in range(num_constraints):
            row = A[i][:]
            slack = [0]*num_constraints
            artificial = [0]*num_constraints
            if signs[i] == "<=":
                # Agregar variable de holgura
                slack[i] = 1
                c_extended.append(0)
                basis.append(var_index)
                var_index += 1
            elif signs[i] == ">=":
                # Agregar variable de superávit y variable artificial
                slack[i] = -1
                artificial[i] = 1
                c_extended.extend([0, M])
                artificial_vars.append(var_index + 1)
                basis.append(var_index + 1)
                var_index += 2
            elif signs[i] == "=":
                # Agregar variable artificial
                artificial[i] = 1
                c_extended.append(M)
                artificial_vars.append(var_index)
                basis.append(var_index)
                var_index += 1
            else:
                self.output_text.insert(tk.END, "Signo de restricción no reconocido.\n")
                return
            row.extend(slack)
            row.extend(artificial)
            row.append(b[i])
            tableau.append(row)
        
        # Ajustar longitud de c_extended
        while len(c_extended) < len(tableau[0])-1:
            c_extended.append(0)
        
        # Si es minimización, cambiar signo de M
        if self.optimization_type.get() == "Minimizar":
            c_extended = [coeff if i not in artificial_vars else coeff - M for i, coeff in enumerate(c_extended)]
        else:
            c_extended = [coeff if i not in artificial_vars else coeff + M for i, coeff in enumerate(c_extended)]
        
        # Mostrar el tableau inicial
        self.output_text.insert(tk.END, "Tableau Inicial (Método de la M Grande):\n")
        self.display_tableau(tableau, c_extended, basis)
        
        # Iniciar iteraciones
        iteration = 0
        max_iterations = 100
        while iteration < max_iterations:
            iteration += 1
            self.output_text.insert(tk.END, f"\nIteración {iteration}:\n")
            
            # Calcular zj y cj - zj
            zj = [0]*(len(tableau[0])-1)
            for i in range(len(basis)):
                for j in range(len(zj)):
                    zj[j] += c_extended[basis[i]] * tableau[i][j]
            cj_zj = [c_extended[j] - zj[j] for j in range(len(zj))]
            
            self.output_text.insert(tk.END, f"Ganancias relativas (cj - zj): {['{0:.2f}'.format(num) for num in cj_zj]}\n")
            
            # Verificar optimalidad
            if all(value <= 1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Se encontró la solución óptima.\n")
                break
            
            # Variable entrante
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] > 1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes. El problema puede tener soluciones múltiples.\n")
                break
            entering = max(entering_candidates, key=lambda j: cj_zj[j])
            self.output_text.insert(tk.END, f"Variable entrante: x{entering+1}\n")
            
            # Calcular razones
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))
            self.output_text.insert(tk.END, f"Razones: {['{0:.2f}'.format(num) if num != float('inf') else 'Inf' for num in ratios]}\n")
            
            # Verificar ilimitado
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "La solución es ilimitada.\n")
                return
            
            # Variable saliente
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: x{basis[leaving]+1}\n")
            
            # Pivoteo
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            # Actualizar base
            basis[leaving] = entering
            
            # Mostrar tableau
            self.display_tableau(tableau, c_extended, basis)
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones.\n")
            return
        
        # Verificar si hay variables artificiales en la base
        if any(var in basis for var in artificial_vars):
            self.output_text.insert(tk.END, "El problema no tiene solución factible.\n")
            return
        
        # Extraer solución
        solution = [0]*(len(tableau[0])-1)
        for i in range(len(basis)):
            solution[basis[i]] = tableau[i][-1]
        
        z = sum(c_extended[i]*solution[i] for i in range(len(c_extended)))
        if self.optimization_type.get() == "Minimizar":
            z = -z
        self.output_text.insert(tk.END, f"\nSolución Óptima:\n")
        for i in range(num_variables):
            self.output_text.insert(tk.END, f"x{i+1} = {solution[i]:.2f}\n")
        self.output_text.insert(tk.END, f"Z = {z:.2f}\n")
        
    def two_phase_method(self, c, A, b, signs, var_constraints):
        self.output_text.delete(1.0, tk.END)
        
        num_variables = len(c)
        num_constraints = len(A)
        
        c_extended = c[:]
        tableau = []
        basis = []
        artificial_vars = []
        var_index = num_variables
        
        # Fase 1: Construir problema auxiliar
        c_phase1 = [0]*(num_variables)
        for i in range(num_constraints):
            row = A[i][:]
            slack = [0]*num_constraints
            artificial = [0]*num_constraints
            if signs[i] == "<=":
                slack[i] = 1
                c_phase1.append(0)
                basis.append(var_index)
                var_index += 1
            elif signs[i] == ">=":
                slack[i] = -1
                artificial[i] = 1
                c_phase1.extend([0, 1])
                artificial_vars.append(var_index + 1)
                basis.append(var_index + 1)
                var_index += 2
            elif signs[i] == "=":
                artificial[i] = 1
                c_phase1.append(1)
                artificial_vars.append(var_index)
                basis.append(var_index)
                var_index += 1
            else:
                self.output_text.insert(tk.END, "Signo de restricción no reconocido.\n")
                return
            row.extend(slack)
            row.extend(artificial)
            row.append(b[i])
            tableau.append(row)
        
        # Ajustar longitud de c_phase1
        while len(c_phase1) < len(tableau[0])-1:
            c_phase1.append(0)
        
        # Mostrar el tableau inicial de la Fase 1
        self.output_text.insert(tk.END, "Fase 1: Tableau Inicial\n")
        self.display_tableau(tableau, c_phase1, basis)
        
        # Iniciar iteraciones de la Fase 1
        iteration = 0
        max_iterations = 100
        while iteration < max_iterations:
            iteration += 1
            self.output_text.insert(tk.END, f"\nFase 1 - Iteración {iteration}:\n")
            
            # Calcular zj y cj - zj
            zj = [0]*(len(tableau[0])-1)
            for i in range(len(basis)):
                for j in range(len(zj)):
                    zj[j] += c_phase1[basis[i]] * tableau[i][j]
            cj_zj = [c_phase1[j] - zj[j] for j in range(len(zj))]
            
            self.output_text.insert(tk.END, f"Ganancias relativas (cj - zj): {['{0:.2f}'.format(num) for num in cj_zj]}\n")
            
            # Verificar optimalidad
            if all(value >= -1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Fase 1 completada.\n")
                break
            
            # Variable entrante
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] < -1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes en Fase 1.\n")
                break
            entering = min(entering_candidates, key=lambda j: cj_zj[j])
            self.output_text.insert(tk.END, f"Variable entrante: x{entering+1}\n")
            
            # Calcular razones
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))
            self.output_text.insert(tk.END, f"Razones: {['{0:.2f}'.format(num) if num != float('inf') else 'Inf' for num in ratios]}\n")
            
            # Verificar ilimitado
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "El problema no tiene solución factible.\n")
                return
            
            # Variable saliente
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: x{basis[leaving]+1}\n")
            
            # Pivoteo
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            # Actualizar base
            basis[leaving] = entering
            
            # Mostrar tableau
            self.display_tableau(tableau, c_phase1, basis)
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones en Fase 1.\n")
            return
        
        # Verificar si valor óptimo es cero
        z0 = sum(c_phase1[basis[i]] * tableau[i][-1] for i in range(len(basis)))
        if abs(z0) > 1e-5:
            self.output_text.insert(tk.END, "El problema no tiene solución factible.\n")
            return
        
        # Eliminar columnas de variables artificiales
        artificial_vars = sorted(artificial_vars, reverse=True)
        for var in artificial_vars:
            for row in tableau:
                del row[var]
            del c_extended[var]
        
        # Actualizar base
        basis = [b for b in basis if b not in artificial_vars]
        
        # Ajustar c_extended
        while len(c_extended) < len(tableau[0])-1:
            c_extended.append(0)
        
        # Fase 2: Resolver el problema original
        self.output_text.insert(tk.END, "\nFase 2: Resolver el problema original\n")
        self.display_tableau(tableau, c_extended, basis)
        
        # Iniciar iteraciones de la Fase 2
        iteration = 0
        while iteration < max_iterations:
            iteration += 1
            self.output_text.insert(tk.END, f"\nFase 2 - Iteración {iteration}:\n")
            
            # Calcular zj y cj - zj
            zj = [0]*(len(tableau[0])-1)
            for i in range(len(basis)):
                for j in range(len(zj)):
                    zj[j] += c_extended[basis[i]] * tableau[i][j]
            cj_zj = [c_extended[j] - zj[j] for j in range(len(zj))]
            
            self.output_text.insert(tk.END, f"Ganancias relativas (cj - zj): {['{0:.2f}'.format(num) for num in cj_zj]}\n")
            
            # Verificar optimalidad
            if all(value <= 1e-5 for value in cj_zj):
                self.output_text.insert(tk.END, "Se encontró la solución óptima.\n")
                break
            
            # Variable entrante
            entering_candidates = [j for j in range(len(cj_zj)) if cj_zj[j] > 1e-5]
            if not entering_candidates:
                self.output_text.insert(tk.END, "No hay variables entrantes. El problema puede tener soluciones múltiples.\n")
                break
            entering = max(entering_candidates, key=lambda j: cj_zj[j])
            self.output_text.insert(tk.END, f"Variable entrante: x{entering+1}\n")
            
            # Calcular razones
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][entering] > 1e-5:
                    ratios.append(tableau[i][-1] / tableau[i][entering])
                else:
                    ratios.append(float('inf'))
            self.output_text.insert(tk.END, f"Razones: {['{0:.2f}'.format(num) if num != float('inf') else 'Inf' for num in ratios]}\n")
            
            # Verificar ilimitado
            if all(r == float('inf') for r in ratios):
                self.output_text.insert(tk.END, "La solución es ilimitada.\n")
                return
            
            # Variable saliente
            leaving = ratios.index(min(ratios))
            self.output_text.insert(tk.END, f"Variable saliente: x{basis[leaving]+1}\n")
            
            # Pivoteo
            pivot_element = tableau[leaving][entering]
            tableau[leaving] = [x / pivot_element for x in tableau[leaving]]
            
            for i in range(len(tableau)):
                if i != leaving:
                    factor = tableau[i][entering]
                    tableau[i] = [tableau[i][j] - factor * tableau[leaving][j] for j in range(len(tableau[i]))]
            
            # Actualizar base
            basis[leaving] = entering
            
            # Mostrar tableau
            self.display_tableau(tableau, c_extended, basis)
        else:
            self.output_text.insert(tk.END, "Se alcanzó el número máximo de iteraciones en Fase 2.\n")
            return
        
        # Extraer solución
        solution = [0]*(len(tableau[0])-1)
        for i in range(len(basis)):
            solution[basis[i]] = tableau[i][-1]
        
        z = sum(c_extended[i]*solution[i] for i in range(len(c_extended)))
        if self.optimization_type.get() == "Minimizar":
            z = -z
        self.output_text.insert(tk.END, f"\nSolución Óptima:\n")
        for i in range(num_variables):
            self.output_text.insert(tk.END, f"x{i+1} = {solution[i]:.2f}\n")
        self.output_text.insert(tk.END, f"Z = {z:.2f}\n")

if __name__ == "__main__":
    root = tk.Tk()
    app = LinearProgrammingSolver(root)
    root.mainloop()
