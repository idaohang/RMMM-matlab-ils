#include<stdio.h>
#include<malloc.h>
void insert();
void del();
void display();

struct node
{
    int priority;
    int info;
    struct node *next;
}*start,*q,*temp,*new;

typedef struct node *N;

int main()
{
    //Just testing the pq
    insert(7,3);
    insert(3,1);
    insert(4,2);
    insert(5,1);
    insert(6,1);
    
    printf("%d\n",pop());
    printf("%d\n",pop());
    printf("%d\n",pop());
    printf("%d\n",pop());
    printf("%d\n",pop());
}

void insert(int item,int itprio)
{
    new=(N*)malloc(sizeof(N));
    new->info=item;
    new->priority=itprio;
    if(start==NULL || itprio <= start->priority)
    {
        new->next=start;
        start=new;
    }
    else
    {
        q=start;
        while(q->next != NULL && q->next->priority<=itprio)
        q=q->next;
        new->next=q->next;
        q->next=new;
    }
}

int pop()
{
    int top;
    if(start==NULL)
    {
        printf("\nQUEUE UNDERFLOW\n");
    }
    else
    {
        new=start;
        top = new->info;
        start=start->next;
        free(start);
    }
    return top;
}

int peek()
{
    return start->info;
} 
