from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button

class DNASequenceComplementApp(App):
    def build(self):
        layout = BoxLayout(orientation='vertical', padding=20, spacing=10)

        label = Label(text="Enter DNA sequence:", font_size=18)
        layout.add_widget(label)

        self.entry = TextInput(font_size=16, size_hint_y=None, height=40)
        layout.add_widget(self.entry)

        calculate_button = Button(text="Calculate", on_press=self.calculate_complements, size_hint_y=None, height=40)
        layout.add_widget(calculate_button)

        self.result_text = TextInput(font_size=16, size_hint_y=None, height=200, readonly=True)
        layout.add_widget(self.result_text)

        return layout

    def complement_base(self, base):
        if base == 'A':
            return 'T'
        elif base == 'T':
            return 'A'
        elif base == 'C':
            return 'G'
        elif base == 'G':
            return 'C'
        else:
            return base  # If the base is not A, T, C, or G, return the same base

    def reverse_complement(self, seq):
        complement_seq = [self.complement_base(base) for base in seq]
        reverse_complement_seq = ''.join(complement_seq[::-1])
        return reverse_complement_seq

    def calculate_complements(self, instance):
        dna_sequence = self.entry.text.upper()

        complement_sequence = ''.join([self.complement_base(base) for base in dna_sequence])
        reverse_complement_sequence = self.reverse_complement(dna_sequence)

        result_text = f"Complement: {complement_sequence}\nReverse complement: {reverse_complement_sequence}"
        self.result_text.text = result_text

if __name__ == '__main__':
    app = DNASequenceComplementApp()
    app.title = "DNA Sequence Complement"   # Set the app's title with spaces
    app.run()
