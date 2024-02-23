using DataStructures

# Graham scan for a set of complex numbers
function Graham_scan(points)
        if length(points) < 4
                return points
        end
        sort!(points, by=imag, rev=true)
        points[2:end] = sort!(points[2:end], by = z -> angle(z - points[1]))

        s = Stack{Complex}()
        for point in points
                while length(s) > 1 && ccw(next_to_top(s), first(s), point) <= 0
                        pop!(s)
                end
                push!(s, point)
        end
        return s
end

# get the second-from-last element from a stack
function next_to_top(s::Stack)
        temp = pop!(s)
        val = first(s)
        push!(s, temp)
        return val
end

# check the orientation of three complex numbers 
function ccw(x::Complex, y::Complex, z::Complex)
        if angle(y - x) < angle(z - x)
                return 1
        elseif angle(y - x) > angle(z - x)
                return -1
        else
                return 0
        end
end

function test()
        points = [0 + im, 1 + im, 1 + 0*im, 0.75 + 0.75*im]
        display(Graham_scan(points))
end
