package main

import (
	"bufio"
	"container/heap"
	"fmt"
	"os"
	"strings"
	"sync"
)

// Структура для хранения профиля белка с биграммами
type ProteinProfile struct {
	Profile [676]uint32
	Sum     uint64
}

// Структура для хранения результатов родственности
type ProteinSimilarity struct {
	similarity float64
	name       string
}

// Приоритетная очередь для топ-100 родственников
type MaxHeap []ProteinSimilarity

func (h MaxHeap) Len() int           { return len(h) }
func (h MaxHeap) Less(i, j int) bool { return h[i].similarity > h[j].similarity }
func (h MaxHeap) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }

func (h *MaxHeap) Push(x interface{}) {
	*h = append(*h, x.(ProteinSimilarity))
}

func (h *MaxHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

// Построение профиля для одной последовательности с биграммами
func buildProfile(sequence string) ProteinProfile {
	var profile [676]uint32
	var sum uint64

	for i := 0; i+1 < len(sequence); i++ {
		index := (int(sequence[i]-'A') * 26) + int(sequence[i+1]-'A')
		profile[index]++
		sum++
	}

	return ProteinProfile{profile, sum}
}

// Вычисление родственности между двумя профилями
func calculateSimilarity(X, Y ProteinProfile) float64 {
	sumMin := uint64(0)
	for i := 0; i < 676; i++ {
		if X.Profile[i] < Y.Profile[i] {
			sumMin += uint64(X.Profile[i])
		} else {
			sumMin += uint64(Y.Profile[i])
		}
	}
	sumMax := X.Sum + Y.Sum - sumMin
	return float64(sumMin) / float64(sumMax)
}

// Построение профилей параллельно с ограничением горутин
func buildProfilesParallel(sequences []string, maxGoroutines int, startIndex int) []ProteinProfile {
	profiles := make([]ProteinProfile, len(sequences))
	var wg sync.WaitGroup
	semaphore := make(chan struct{}, maxGoroutines)

	for i, seq := range sequences {
		wg.Add(1)
		semaphore <- struct{}{}
		go func(index int, sequence string) {
			defer wg.Done()
			profile := buildProfile(sequence)
			profiles[index] = profile
			<-semaphore
		}(i, seq)
	}

	wg.Wait()
	return profiles
}

// Поиск топ-100 родственников с логированием индексов
func findTop100Similar(newProteinSequence string, profiles []ProteinProfile, names []string, topN int) []ProteinSimilarity {
	newProfile := buildProfile(newProteinSequence)
	var minHeap MaxHeap
	heap.Init(&minHeap)

	for i, profile := range profiles {
		similarity := calculateSimilarity(newProfile, profile)
		heap.Push(&minHeap, ProteinSimilarity{-similarity, names[i]})
		if minHeap.Len() > topN {
			heap.Pop(&minHeap)
		}
	}

	var result []ProteinSimilarity
	for minHeap.Len() > 0 {
		item := heap.Pop(&minHeap).(ProteinSimilarity)
		item.similarity = -item.similarity
		result = append(result, item)
	}

	return result
}

// Обработка файла по частям
func processFileInChunks(filename string, chunkSize int, maxGoroutines int) ([]string, []ProteinProfile, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, fmt.Errorf("ошибка открытия файла: %v", err)
	}
	defer file.Close()

	var allNames []string
	var allProfiles []ProteinProfile
	scanner := bufio.NewScanner(file)
	var currentChunkNames []string
	var currentChunkSequences []string
	var currentName string
	var currentSequence strings.Builder
	var lineCount int
	var totalProcessed int // Общее количество обработанных белков

	for scanner.Scan() {
		lineCount++
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if currentName != "" {
				currentChunkNames = append(currentChunkNames, currentName)
				currentChunkSequences = append(currentChunkSequences, currentSequence.String())
				currentSequence.Reset()
			}
			currentName = line
		} else {
			currentSequence.WriteString(line)
		}

		if len(currentChunkSequences) >= chunkSize {
			profiles := buildProfilesParallel(currentChunkSequences, maxGoroutines, totalProcessed)
			allNames = append(allNames, currentChunkNames...)
			allProfiles = append(allProfiles, profiles...)
			totalProcessed += len(currentChunkSequences)
			currentChunkNames = nil
			currentChunkSequences = nil
		}
	}

	if currentName != "" {
		currentChunkNames = append(currentChunkNames, currentName)
		currentChunkSequences = append(currentChunkSequences, currentSequence.String())
	}

	if len(currentChunkSequences) > 0 {
		profiles := buildProfilesParallel(currentChunkSequences, maxGoroutines, totalProcessed)
		allNames = append(allNames, currentChunkNames...)
		allProfiles = append(allProfiles, profiles...)
	}

	return allNames, allProfiles, scanner.Err()
}

func main() {
	// Создаём файл для вывода результатов
	outputFile, err := os.Create("output.txt")
	if err != nil {
		fmt.Printf("Не удалось создать файл вывода: %v\n", err)
		fmt.Println("Нажмите Enter для выхода...")
		fmt.Scanln()
		return
	}
	defer outputFile.Close()

	// Запрашиваем у пользователя название файла
	var filename string
	fmt.Print("Введите название файла с данными: ")
	fmt.Scanln(&filename)

	// Обрабатываем файл по частям, с ограничением горутин
	names, profiles, err := processFileInChunks(filename, 10000, 6)

	fmt.Printf("Успешно прочитано %d белков.\n", len(names))
	outputFile.WriteString(fmt.Sprintf("Успешно прочитано %d белков.\n", len(names)))

	// Запрашиваем у пользователя последовательность белка
	var newProteinSequence string
	fmt.Print("Введите последовательность белка для сравнения: ")
	scanner := bufio.NewScanner(os.Stdin)
	if scanner.Scan() {
		newProteinSequence = scanner.Text()
	}

	//newProteinSequence := "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL"
	//newProteinSequence := "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"

	topN := 100
	topSimilarities := findTop100Similar(newProteinSequence, profiles, names, topN)

	// Выводим топ-100 родственников в файл
	outputFile.WriteString(fmt.Sprintf("Toп %d родственников:\n", topN))
	for i := len(topSimilarities) - 1; i >= 0; i-- {
		line := fmt.Sprintf("%d. %s (родственность: %.4f)\n", len(topSimilarities)-i, topSimilarities[i].name, topSimilarities[i].similarity)
		outputFile.WriteString(line)
		fmt.Print(line)
	}

	fmt.Println("Результаты записаны в файл output.txt")
}
