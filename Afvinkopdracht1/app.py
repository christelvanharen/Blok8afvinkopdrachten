from flask import Flask, render_template, request
import mysql.connector

app = Flask(__name__)


@app.route('/', methods=["POST", "GET"])
def webpagina():
    """
    Zorgt ervoor dat de webapplicatie wordt weergeven dmv de html
    :return: de webapplicatie
    """
    if request.method == "POST":
        zoek = request.form.get("zoek", "")
        rows = connect_database(zoek)

        return render_template("afvink1.html", database=rows, zoek=zoek)

    else:
        rows = connect_database("None")
        return render_template("afvink1.html", database=rows,
                               zoek="None")


def connect_database(zoek):
    """
    Haalt de description uit de ensembldb database
     en filtert deze op het zoekwoord
    :param zoek: Het ingegeven zoekwoord
    :return: Een lijst met de juiste discriptions
    """
    print("Het huidige zoekwoord is:", zoek)
    conn = mysql.connector.connect(host='ensembldb.ensembl.org',
                                   user='anonymous',
                                   db='homo_sapiens_core_95_38')
    cursor = conn.cursor()
    cursor.execute("select description from gene")
    rows = cursor.fetchall()
    des = []
    for row in rows:
        if str(row) != "(None,)":
            if zoek.upper() in str(row).upper():

                regel = []
                door = True
                while door:
                    pos_start = str(row).upper().index(zoek.upper())
                    pos_end = pos_start + len(zoek)
                    woord = str(row)[pos_start:pos_end]
                    regel.append(str(row)[0:pos_start])
                    regel.append(woord)
                    row = str(row)[pos_end:]
                    door = zoek.upper() in str(row).upper()

                regel.append(row)
                des.append(regel)

            elif zoek == 'None':
                des.append(row)

    cursor.close()
    conn.close()
    return des


if __name__ == '__main__':
    app.run()