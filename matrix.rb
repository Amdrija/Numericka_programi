class Matrix
  attr_reader :matrix
  def initialize(dimension)
    @dimension = dimension
    @matrix = []
    dimension.times do |i|
      @matrix.push([])
      dimension.times do |j|
        @matrix[i][j] = 0
      end
    end
  end

  def read
    @dimension.times do |i|
      @dimension.times do |j|
        @matrix[i][j] = gets.to_f
      end
    end
  end

  def read_from_file(file)
    file_data = file.readlines
    @dimension.times do |i|
      @dimension.times do |j|
        @matrix[i][j] = file_data[j + i * @dimension].to_f
      end
    end
    file_data
  end

  def to_s(precision = 4)
    @dimension.times do |i|
      @dimension.times do |j|
        print @matrix[i][j].round(precision)
        print ", " if j != @dimension - 1
      end
      print "\n"
    end
  end

  def swap(i, j, option_hash = {row: true})
    #TODO: add checking if i or j are in dimensions
    if(option_hash[:row])
      swap_row(i, j)
    else
      swap_column(i, j)
    end
  end

  def swap_row(i, j)
    #TODO: add checking if i and j are in dimensions
    helper = @matrix[i]
    @matrix[i] = @matrix[j]
    @matrix[j] = helper
  end

  def swap_column(i, j)
    #TODO: add checking if i and j are in dimensions
    @dimension.times do |row|
      helper = @matrix[row][i]
      @matrix[row][i] = @matrix[row][j]
      @matrix[row][j] = helper
    end
  end

  def add_row_to_row(i, j, coef = 1)
    #TODO: add checking if i and j are in dimensions
    @dimension.times do |column|
      @matrix[j][column] += coef * @matrix[i][column]
    end
  end

  def add_column_to_column(i, j, coef = 1)
    #TODO: add checking if i and j are in dimensions
    @dimension.times do |row|
      @matrix[row][j] += coef * @matrix[row][i]
    end
  end

  def swap_row_with_maximum(i)
    maximum_row = find_maximum(i)
    swap_row(i, maximum_row) if i != maximum_row
  end

  def find_maximum(n)
    maximum_row_index = n
    maximum_row_abs_value = @matrix[n][n].abs
    (@dimension - n - 1).times do |i|
      #the row is j + i + 1, because we need the next row to i-th
      #and that's i + 1, i + 2 etc.
      if @matrix[n + i + 1][n].abs > maximum_row_abs_value 
        maximum_row_abs_value = @matrix[n + i + 1][n]
        maximum_row_index = n + i + 1
      end
    end
    maximum_row_index
  end
end

class LinearSystemOfEquations
  def initialize(dimension, precision)
    @matrix = Matrix.new(dimension)
    @dimension = dimension
    @free_coef = []
    @dimension.times do |i|
      @free_coef.push(0)
    end
    @precision = precision
  end

  def read
    read_matrix
    read_free_coef
  end

  def read_from_file(file)
    file_data = @matrix.read_from_file(file)

    @dimension.times do |i|
      @free_coef[i] = file_data[@dimension * @dimension + i].to_f
    end
  end

  def read_matrix
    @matrix.read
  end

  def read_free_coef
    @dimension.times do |i|
      @free_coef[i] = gets.chomp.to_f
    end
  end

  def to_s
    @dimension.times do |i|
      @dimension.times do |j|
        print @matrix.matrix[i][j].round(@precision)
        print ", " if j != @dimension - 1
      end
      print " = "
      print @free_coef[i].round(@precision)
      print "\n"
    end
  end

  def gaus_solve(pivot = true)
    lowwer_triangle_matrix(pivot)
    reverse_solve
  end

  def lowwer_triangle_matrix(pivot_flag = true)
    @dimension.times do |i|
      #TODO: Add checking if the matrix[i][i] == 0
      all_rows_0 = false #true if the i-th coeficient of all rows are 0
      #if it's true, that mens we don't have to do gaus_method
      if pivot_flag
        maximum_row = @matrix.find_maximum(i)
        if maximum_row != i
          @matrix.swap_row(i, maximum_row)
          helper = @free_coef[i]
          @free_coef[i] = @free_coef[maximum_row]
          @free_coef[maximum_row] = helper
        end
      end
      (@dimension - i - 1).times do |j|
        coef = -  @matrix.matrix[j + i + 1][i] * 1.0 / @matrix.matrix[i][i]
        @matrix.add_row_to_row(i, j + i + 1, coef)
        #adds i-th row multiplied by the caslculated coefficient to (j + i + 1)-th row 
        #the other rows are offset with j + i + 1, because we the inner loop is starting
        #from 0, but we want it to start from next row to i-th , and that's (i + 1)-th row
        @free_coef[j + i + 1] += coef * @free_coef[i]
      end
    end
    puts self
  end

  def reverse_solve
    solution = []
    (@dimension).times do |i|
      position = @dimension - 1 - i
      solution.unshift((@free_coef[position] - sum_up_to_n(position, solution)) / @matrix.matrix[position][position])
      #according to formula, we add the free coefficient at that position with
      #all the other solution variables ,that are known up until then, multiplied
      #by the coefficient at their positoon in the row
    end
    solution.map {|el| el.round(@precision)}
  end

  def check_solution(solution)
    @matrix.matrix.each_with_index do |row, i|
      left_side = 0
      row.each_with_index do |coef, j|
        left_side += coef * solution[j]
      end
      puts left_side - @free_coef[i]
    end
  end

  private
  #current solution is incomplete
  #it contains all known variables from
  #system dimension to n
  def sum_up_to_n(n, current_solution)
    sum = 0
    current_solution.each_with_index do |variable, i|
      sum += variable * @matrix.matrix[n][i + n + 1]
      #the column index is i + n + 1, because we need
      #to offset the current solution index with n + 1
      #missing solution variables
    end
    sum
  end

end

dimension = 4
precision = 6
sys_pivot = LinearSystemOfEquations.new(dimension, precision)
sys_no_pivot = LinearSystemOfEquations.new(dimension, precision)
sys_original = LinearSystemOfEquations.new(dimension, precision)
file = File.open("test_domaci8.txt")
sys_pivot.read_from_file(file)
file.rewind
sys_no_pivot.read_from_file(file)
file.rewind
sys_original.read_from_file(file)
file.close
puts sys_original
solution_pivot = sys_pivot.gaus_solve(true)
solution_no_pivot = sys_no_pivot.gaus_solve(false)
puts solution_pivot
puts solution_no_pivot
sys_pivot.check_solution(solution_pivot)
sys_no_pivot.check_solution(solution_no_pivot)
